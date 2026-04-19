library(tidyverse)
library(tools)
library(terra)
setwd('D:/Rwork/connectivity/aqua_data')

process_species_data <- function(species_name) {
  depth_layers <- list(
    c(0,5),c(5,10),c(10,25),c(25,50),c(50,100),c(100,150),c(150,200)
  )
  raw_dir <- file.path('.', species_name, paste0(species_name, '_raw'))
  aqua_files <- list.files(path = raw_dir,
                           pattern = '\\.csv$',
                           full.names = TRUE,
                           recursive = TRUE)
  
  for(d in seq_along(depth_layers)){
    range_i <- depth_layers[[d]]
    dmin <- range_i[1]
    dmax <- range_i[2]
    depth_label <- paste(dmin,dmax,sep = '_')
    result_df <- data.frame(
      Genus = character(0),
      Species = character(0),
      Latitude = numeric(0),
      Longitude = numeric(0),
      depth = character(0),
      stringsAsFactors = FALSE
    )
    for(f in aqua_files){
      data <- read.csv(f,skip = 8,header = FALSE,stringsAsFactors = FALSE)
      depth_val <- as.numeric(data[, 6])
      depth_index <- which(depth_val >= dmin & depth_val < dmax)
      if(length(depth_index) == 0) next
      depth_layer_data <- data[depth_index,]
      tmp_df <- data.frame(
        Genus = depth_layer_data[,1],
        Species = depth_layer_data[,2],
        Latitude = as.numeric(depth_layer_data[,3]),
        Longitude = as.numeric(depth_layer_data[,4]),
        depth = depth_label,
        stringsAsFactors = FALSE
      )
      result_df <- rbind(result_df, tmp_df)
    }
    result_df$occurance <- rep(1,nrow(result_df))
    clear_dir <- file.path('.', species_name, species_name)
    if (!dir.exists(clear_dir)) {
      dir.create(clear_dir, recursive = TRUE)
    }
    output_name <- file.path(clear_dir, 
                             paste0(species_name, "_data_", depth_label, "m_2003.csv"))
    write.csv(result_df, file = output_name, row.names = FALSE)
  }
  
  # rasterize
  rasterize_func <- function(dt,value){
    v <- vect(dt,geom = c('Longitude','Latitude'),crs = 'EPSG:4326')
    r_template <- rast(res = 0.5,
                       ext = ext(-180,180,-90,90),
                       crs = 'EPSG:4326')
    rs <- rasterize(v,r_template,field = value,fun = 'first')
    return(rs)
  }
  
  clear_dir <- file.path('.', species_name, species_name)
  read_pattern <- paste0("^", species_name, "_data_.*\\.csv$")
  clear_aqua_files <- list.files(path = clear_dir,
                                 pattern = read_pattern,
                                 full.names = TRUE,
                                 recursive = FALSE)
  out_dir <- file.path('.', species_name, paste0(species_name, "_rasterize"))
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  for(f in clear_aqua_files){
    dt <- read.csv(f)
    fname <- basename(f)
    fname_rast <- sub('\\.csv$','.tif',fname)
    rs <- rasterize_func(dt,'occurance')
    out_path <- file.path(out_dir, fname_rast)
    writeRaster(rs,out_path,overwrite = TRUE)
    
  }
}

process_species_data('Megaptera')
process_species_data('Hippoglossus')
process_species_data('Gadus')

# 重采样见Python   


