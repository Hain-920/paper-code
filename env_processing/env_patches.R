library(hdf5r)
library(R.matlab)
library(raster)
library(sf)
library(sp)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(curl)
library(terra)
library(stringr)
library(lwgeom)
library(devtools)
library(tidyterra)
library(purrr)
library(viridis)
library(stringr)
library(ncdf4)
library(magrittr)

# 溶解氧导入
pointer <- H5File$new("do_predict_5902m_softmax_prior_rmse.mat", mode = "r")
pointer$ls(recursive = TRUE)
lon <- pointer[["/lon"]]$read()
lat <- pointer[["/lat"]]$read()
depth <- pointer[['/depth']]$read()
y_data <- pointer[['/y']]$read()
pointer$close_all()
dim(y_data)

month_index = 1:12
time_index = 44

# 温度导入
data_dir <- 'T_2003'
t_files <- list.files(
  data_dir,
  pattern = '^T_regrid_2003[0-9]{2}\\.mat$',
  full.names = TRUE
)

for(i in seq_along(t_files)){
  f <- H5File$new(t_files[i], mode = "r")
  y <- f[['/T_regrid']]$read()
  f$close_all()
  if(i==1){
    T_sum <- ifelse(is.na(y),0,y)
    n <- !is.na(y)
  }else{
    T_sum <- T_sum + ifelse(is.na(y),0,y)
    n <- n + !is.na(y)
  }
  
}
T_2003_mean <- T_sum/n
lon <- seq(0.5,359.5,by=1)
lat <- seq(-89.5,89.5,by=1)
writeMat('T_2003/T_2003_mean.mat',
         T_2003 = T_2003_mean,
         lon = lon,
         lat = lat,
         depth = depth
)

# 盐度导入
sal_data_dir <- 'S_regrid'
sal_files <- list.files(
  sal_data_dir,
  pattern = '^S_regrid_2003[0-9]{2}\\.mat$',
  full.names = TRUE
)

for(i in seq_along(sal_files)){
  f <- H5File$new(sal_files[i], mode = "r")
  y <- f[['/S_regrid']]$read()
  f$close_all()
  if(i==1){
    S_sum <- ifelse(is.na(y),0,y)
    n <- !is.na(y)
  }else{
    S_sum <- S_sum + ifelse(is.na(y),0,y)
    n <- n + !is.na(y)
  }
  
}
sal_2003_mean <- S_sum/n
lon <- seq(0.5,359.5,by=1)
lat <- seq(-89.5,89.5,by=1)
writeMat('sal_2003/sal_2003_mean.mat',
         sal_2003 = sal_2003_mean,
         lon = lon,
         lat = lat,
         depth = depth
)

# pH导入
ph <- rast('pH_data/DATA2023000009_2003.nc',subds='pH')
ph[ph==-999] <- NA
ph_year_depth <- tapp(
  ph,   # 数据
  index = rep(1:41,each=12),   # 索引
  fun = mean,
  na.rm = TRUE
)
ph_0_200 <- ph_year_depth[[1:17]]
depths <- c(0,5,10,20,30,40,50,60,70,80,90,100,120,140,160,180,200)
for(i in 1:17){
  ph_depth <- ph_0_200[[i]]
  depth_idx <- depths[i]
  outname <- sprintf('pH_data/pH_2003_%dm.tif',depth_idx)
  writeRaster(ph_depth,outname, overwrite = TRUE)
}
writeRaster(ph_year_depth,'pH_data/pH_2003.tif',overwrite = TRUE)


# 深度分层输出文件
depth_layers <- list(
  c(0,5),c(5,10),c(10,25),c(25,50),c(50,100),c(100,150),c(150,200)
)
# DO
for(i in seq_along(depth_layers)){
  range_i <- depth_layers[[i]]
  dmin <- range_i[1]
  dmax <- range_i[2]
  depth_index <- which(depth>=dmin & depth<dmax)
  DO_layer <- apply(y_data[,,depth_index,month_index,time_index,drop=FALSE],c(1,2),mean,na.rm=TRUE)
  DO_layer <- t(DO_layer)
  DO_layer <- DO_layer[nrow(DO_layer):1,]
  DO_layer <- rast(DO_layer) 
  DO_layer <- rotate(DO_layer)
  ext(DO_layer) <- c(min(lon),max(lon),min(lat),max(lat))
  crs(DO_layer) <- '+proj=longlat +datum=WGS84'
  first_col <- DO_layer[,1]
  last_col <- DO_layer[,ncol(DO_layer)]
  avg_col <- (first_col+last_col)/2
  DO_layer[,ncol(DO_layer)] <- avg_col
  output_tif <- sprintf("%s_%d_%dm_%d.tif", 'DO_layer', dmin, dmax,1959+time_index) 
  writeRaster(DO_layer,output_tif,overwrite = TRUE)
  
}

# 温度
data <- readMat("T_2003/T_2003_mean.mat")
t_data <- data$T.2003
t_lon <- data$lon
t_lat <- data$lat
t_depth <- data$depth
for(i in seq_along(depth_layers)){
  range_i <- depth_layers[[i]]
  dmin <- range_i[1]
  dmax <- range_i[2]
  depth_index <- which(t_depth>=dmin & t_depth<dmax)
  t_layer <- apply(t_data[,,depth_index,drop=FALSE],c(1,2),mean,na.rm=TRUE)
  t_layer <- t(t_layer)
  t_layer <- t_layer[nrow(t_layer):1,]
  t_layer <- rast(t_layer) 
  t_layer <- rotate(t_layer)
  ext(t_layer) <- c(-180,180,min(t_lat),max(t_lat))
  crs(t_layer) <- '+proj=longlat +datum=WGS84'
  first_col <- t_layer[,1]
  last_col <- t_layer[,ncol(t_layer)]
  avg_col <- (first_col+last_col)/2
  t_layer[,ncol(t_layer)] <- avg_col
  output_tif <- sprintf("T_layers/%s_%d_%dm_%d.tif", 't_layer', dmin, dmax,2003) 
  writeRaster(t_layer,output_tif,overwrite = TRUE)
}

# 盐度
sal_dataset <- readMat("sal_2003/sal_2003_mean.mat")
sal_data <- sal_dataset$sal.2003
sal_lon <- sal_dataset$lon
sal_lat <- sal_dataset$lat
sal_depth <- sal_dataset$depth
for(i in seq_along(depth_layers)){
  range_i <- depth_layers[[i]]
  dmin <- range_i[1]
  dmax <- range_i[2]
  depth_index <- which(sal_depth>=dmin & sal_depth<dmax)
  sal_layer <- apply(sal_data[,,depth_index,drop=FALSE],c(1,2),mean,na.rm=TRUE)
  sal_layer <- t(sal_layer)
  sal_layer <- sal_layer[nrow(sal_layer):1,]
  sal_layer <- rast(sal_layer) 
  sal_layer <- rotate(sal_layer)
  ext(sal_layer) <- c(-180,180,min(sal_lat),max(sal_lat))
  crs(sal_layer) <- '+proj=longlat +datum=WGS84'
  first_col <- sal_layer[,1]
  last_col <- sal_layer[,ncol(sal_layer)]
  avg_col <- (first_col+last_col)/2
  sal_layer[,ncol(sal_layer)] <- avg_col
  output_tif <- sprintf("sal_layers/%s_%d_%dm_%d.tif", 'sal_layer', dmin, dmax,2003) 
  writeRaster(sal_layer,output_tif,overwrite = TRUE)
}
template <- rast('T_layers/t_layer_0_5m_2003.tif')
sal_dir <- list.files(
  'sal_layers',
  pattern = '_2003\\.tif$',
  full.names = TRUE
)
for(s in sal_dir){
  sal_layer <- rast(s)
  sal_layer_re <- terra::resample(sal_layer,template,method = 'bilinear')
  writeRaster(sal_layer_re,s,overwrite = TRUE)
  
}

# pH
pH_files <- list.files(
  'pH_data',
  pattern = "\\.tif$",
  full.names = TRUE
)
depth <- as.numeric(
  sub(".*_(\\d+)m\\.tif$", "\\1", basename(pH_files))
)   # 提取对应深度
for(i in seq_along(depth_layers)){
  range_i <- depth_layers[[i]]
  dmin <- range_i[1]
  dmax <- range_i[2]
  depth_index <- which(depth>=dmin & depth<dmax)
  ph_layers <- rast(pH_files[depth_index])
  ph_mean <- mean(ph_layers,na.rm = TRUE)
  ph <- t(ph_mean)
  ph <- flip(ph,direction='vertical')
  ph <- rotate(ph)
  crs(ph) <- '+proj=longlat +datum=WGS84'
  output_tif <- sprintf("pH_layers/%s_%d_%dm_%d.tif", 'ph_layer', dmin, dmax,2003) 
  writeRaster(ph,output_tif,overwrite = TRUE)
}

