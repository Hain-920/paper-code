library(terra)
library(magrittr)
library(sf)
setwd('D:/Rwork/connectivity/aqua_data')

# 计算物种在当前环境的cc
climate_change_cal <- function(data,para,lim_lower,lim_upper){
  data$cc <- 0
  for(i in 1:nrow(data)){
    genus_name <- data[i,'Genus']
    spe_name <- data[i,'Species']   # 获取每一行物种名称
    current_lower <- lim_lower[lim_lower[,'Species']==spe_name &
                                 lim_lower[,'Genus']==genus_name,para]   
    current_upper <- lim_upper[lim_upper[,'Species']==spe_name & 
                                 lim_upper[,'Genus']==genus_name,para]   # 获取对应物种的相关环境耐受值
    val <- data[i,para]   # 当前物种所在位置的对应环境值
    if(val < current_lower){
      data$cc[i] <- (current_lower - val)/(current_upper - current_lower)
    }
    if(val > current_upper){
      data$cc[i] <- (val - current_upper)/(current_upper - current_lower)
    }
  }
  return(data)
}

calculate_species_ci <- function(species_name,distance) {
  env_names <- c('DO','pH','Sal','T')
  depth_layers <- c('0_5m','5_10m','10_25m','25_50m','50_100m','100_150m','150_200m')
  weight_mean = 1/4
  lower_sd <- read.csv('lower_sd.csv')
  upper_sd <- read.csv('upper_sd.csv')
  
  out_dir <- file.path('.', species_name, paste0(species_name, '_ci'))
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  for(d in depth_layers){
    csv_path <- file.path('.', species_name, species_name, paste0(species_name, '_data_', d, '_2003.csv'))
    if (!file.exists(csv_path)) {
      next
    }
    aqua_data <- read.csv(csv_path)
    tif_path <- file.path('.', species_name, paste0(species_name, '_rescale'), paste0(species_name, '_data_', d, '_2003.tif'))
    aqua_rasterize <- rast(tif_path)
    aqua_pts <- vect(aqua_data,geom = c("Longitude", "Latitude"),crs = 'EPSG:4326')
    env_stack <- rast(c(
      paste0('env/DO/DO_layer_',d,'_2003.tif'),
      paste0('env/pH/ph_layer_',d,'_2003.tif'),
      paste0('env/sal/sal_layer_',d,'_2003.tif'),
      paste0('env/T/t_layer_',d,'_2003.tif')
    ))   # 环境堆叠栅格
    names(env_stack) <- env_names
    template <- env_stack[[1]]
    R_km <- distance
    R_m <- R_km *1000
    template_proj <- project(template,"EPSG:6933")
    aqua_pts_proj <- project(aqua_pts,crs(template_proj))
    dist_rast <- terra::distance(template_proj,aqua_pts_proj)
    accessible_mask <- dist_rast <= R_m   # 筛选小于最大迁移距离的格点
    accessible_mask[accessible_mask != 1] <- NA
    accessible_mask <- project(accessible_mask,env_stack)
    if(sum(!is.na(values(accessible_mask))) == 0){
      env_accessible <- mask(env_stack,aqua_rasterize)
    }else{env_accessible <- mask(env_stack,accessible_mask)}
    pts <- as.points(env_accessible[[1]],values = FALSE)
    pts_xy <- crds(pts)
    nearest_idx <- terra::nearest(pts,aqua_pts)
    nearest_id <- values(nearest_idx)$to_id
    species_id <- aqua_pts$Species[nearest_id]
    genus_id <- aqua_pts$Genus[nearest_id]
    env_values <- terra::extract(env_stack,pts_xy)
    dt <- data.frame(
      Longitude = pts_xy[,1],
      Latitude = pts_xy[,2],
      Genus = genus_id,
      Species = species_id,
      env_values
    )
    dt <- dt[complete.cases(dt[, c("DO", "pH",'Sal',"T")]), ]
    if (nrow(dt) == 0) {
      message(paste("Depth layer", d, "has no valid data, skipping..."))
      next
    }
    lonlat <- dt[,c('Longitude','Latitude')]
    for(p in env_names){
      cal <- climate_change_cal(dt,p,lower_sd,upper_sd)
      lonlat[[paste0("cc_", p)]] <- cal$cc
    }   
    colnames(lonlat)[3:6] <- env_names
    lonlat$CI_total <- 
      lonlat$DO * weight_mean +
      lonlat$pH * weight_mean +
      lonlat$Sal * weight_mean +
      lonlat$T * weight_mean
    v <- vect(lonlat,geom = c('Longitude','Latitude'),crs ='EPSG:4326')
    
    out_r <- rast(
      resolution = 1,
      extent = ext(-180, 180, -90, 90),
      crs = "EPSG:4326"
    )
    out_r <- rasterize(
      v,
      out_r,
      field = 'CI_total',
      fun = 'sum'
    )
    out_path <- file.path(out_dir, paste0(species_name, '_ci_', d, '_2003.tif'))
    writeRaster(out_r, out_path, overwrite = TRUE)
  }
}
calculate_species_ci('Megaptera',5000)
calculate_species_ci('Hippoglossus',1000)
calculate_species_ci('Gadus',500)
