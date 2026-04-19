library(tools)
library(tidyterra)
library(stringr)
library(dplyr)
library(ggplot2)
setwd('D:/Rwork/connectivity')

# 分级
classification <- function(tif_file){
  env_rast <- rast(tif_file)
  dt <- as.data.frame(env_rast, xy = TRUE, na.rm = TRUE)
  colnames(dt)[3] <- "val"
  dt1 <- dt %>% 
    filter(val == 0) %>% 
    mutate(rank = 0)
  dt2 <- dt[dt$val != 0, ]
  dt2 <- dt2[order(dt2$val), ]
  n <- nrow(dt2)
  dt2$rank <- 1
  dt2[(n/4 + 1):(n/2), "rank"]       <- 2
  dt2[(n/2 + 1):(n*3/4), "rank"]     <- 3
  dt2[(n*3/4 + 1):n, "rank"]         <- 4
  dt_all <- rbind(dt1, dt2)
  pts <- vect(dt_all,geom = c("x", "y"),crs = "EPSG:4326")
  rf <- rast(extent = ext(-180, 180, -90, 90),resolution = c(1, 1),crs = "EPSG:4326")
  r_class <- rasterize(pts,rf,field = "rank",fun = "first",background = NA)
  outfile <- file.path("class_rast",
                       paste0(tools::file_path_sans_ext(basename(tif_file)),"_class.tif")
  )
  writeRaster(r_class, outfile, overwrite=T)
}

# 根据等级赋予板块id
process_species_patches <- function(species_name, year = 2003) {
  depth_layers <- list(
    c(0,5),c(5,10),c(10,25),c(25,50),c(50,100),c(100,150),c(150,200)
  ) 
  
  base_dir <- file.path('aqua_data', species_name)
  ci_dir <- file.path(base_dir, paste0(species_name, '_ci'))
  class_dir <- file.path('class_rast',species_name)
  
  patch_base <- file.path('patches',species_name)
  
  patch_id_dir <- file.path(patch_base, 'patches_id')
  patch_id2rs_dir <- file.path(patch_base, 'patches_id2rs')
  patch_newid_dir <- file.path(patch_base, 'patches_newid')
  patch_elim_rast_dir <- file.path(patch_base, 'patches_elim', 'rast')
  patch_elim_sf_dir <- file.path(patch_base, 'patches_elim', 'sf')
  
  dirs_to_create <- c(class_dir, patch_id_dir, patch_id2rs_dir, 
                      patch_newid_dir, patch_elim_rast_dir, patch_elim_sf_dir)
  for (d in dirs_to_create) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  }
  
  for(d in depth_layers){
    dmin <- d[1]
    dmax <- d[2]
    r_file <- file.path(ci_dir, sprintf("%s_ci_%d_%dm_%d.tif", species_name, dmin, dmax, year))
    t_file <- sprintf("T_layers/t_layer_%d_%dm_%d.tif",dmin, dmax, year)
    if (!file.exists(r_file)) {
      next
    }
    r_class <- classification(r_file)
    t_class <- classification(t_file)
    chi_class <- classification('chi_resample/chi_2003_re.tif')
    t_class_mask <- mask(t_class,r_class)
    chi_class_mask <- mask(chi_class,r_class)
    patch_attr <- t_class * 100 + r_class * 10 + chi_class * 1
    outfile <- file.path(class_dir, sprintf('patch_%d_%dm_%d.tif', dmin, dmax, year))
    writeRaster(patch_attr, outfile, overwrite = TRUE)
  }

  patch_files <- list.files(class_dir,
                            pattern = "^patch_.*\\.tif$", 
                            recursive = TRUE, 
                            full.names = TRUE)
  
  for (fpath in patch_files) {
    f <- basename(fpath)
    depth_range <- str_extract(f, "(?<=patch_).*?(?=_[0-9]{4})")
    fname <- tools::file_path_sans_ext(f)
    
    r <- rast(fpath)
    p <- as.polygons(r, dissolve = TRUE)
    p <- p[!is.na(p$first), ]
    p_single <- terra::disagg(p)
    p_single$pid <- seq_len(nrow(p_single))   # 为每个多边形分配唯一id
    
    writeVector(p_single,file.path(patch_id_dir, paste0(fname, ".shp")),
                overwrite = TRUE)
    r_out <- rast(extent = ext(-180, 180, -90, 90),resolution = 1,crs = crs(r))
    r_pid <- rasterize(p_single, r_out, field = "pid")
    writeRaster(r_pid,file.path(patch_id2rs_dir, f), overwrite = TRUE)
  }
  
  id_files <- list.files(patch_id2rs_dir, pattern = '.tif$', recursive = TRUE, full.names = TRUE)
  
  # 处理合并边界id
  for (f in id_files){
    basefile <- basename(f)
    f_class <- file.path(class_dir, basefile)
    r_patch <- rast(f)
    r_class <- rast(f_class)
    rs <- c(r_class, r_patch)
    names(rs) <- c("class", "patch")   # 合并class和斑块
    dt_all <- as.data.frame(rs, xy=TRUE)
    if(nrow(dt_all)==0)next
    
    dt_left <- subset(dt_all, x == -179.5)
    dt_right <- subset(dt_all, x == 179.5)
    if (nrow(dt_left) == 0 || nrow(dt_right) == 0) {
      next
    }
    row.names(dt_left) <- seq(1, nrow(dt_left))
    row.names(dt_right) <- seq(1, nrow(dt_right))
    
    id = -1
    for (r in 1:nrow(dt_left)){
      left_class <- dt_left$class[r]
      right_class <- dt_right$class[r]   # 提取左右类别
      
      if (is.na(left_class)==FALSE & is.na(right_class)==FALSE){
        if (left_class==right_class){
          dt_all[which(dt_all$patch==dt_left$patch[r]), 'patch'] <- id
          dt_all[which(dt_all$patch==dt_right$patch[r]), 'patch'] <- id
          id = id - 1
        }   # 左右类别相同,合并成新id
      }
    }
    dt_all <- dt_all[complete.cases(dt_all),]
    r_out <- rast(extent = ext(-180,180,-90,90),resolution = 1, crs = crs(rs))
    r_new <- rasterize(vect(dt_all,geom = c('x','y'),crs = "EPSG:4326"),
                       r_out,
                       field = 'patch',
                       fun = 'first')
    poly <- as.polygons(r_new,dissolve = TRUE)
    poly <- st_as_sf(poly)
    colnames(poly) <- c('patch','geometry')
    area_df <- dt_all %>%
      count(patch, name = "shape_area")
    poly <- left_join(poly,area_df,by = 'patch')
    
    outname <- sub("\\.tif$", ".shp", basefile)
    st_write(poly, file.path(patch_newid_dir, outname), append = FALSE, quiet = TRUE)
  }
  
  # 去除小斑块
  patches_id2rs_files <- list.files(patch_id2rs_dir, pattern ="\\.tif$", full.names = TRUE)

  for (fpath in patches_id2rs_files) {
    f <- basename(fpath)
    fname <- tools::file_path_sans_ext(f)
    
    r_pid <- rast(fpath)
    cell_count <- as.data.frame(freq(r_pid, usena = FALSE))
    cell_count <- cell_count[,-1]
    colnames(cell_count) <- c("pid", "ncell")
    small_ids <- cell_count$pid[cell_count$ncell < 9]
    r_clean <- r_pid
    r_clean[r_clean %in% small_ids] <- NA
    
    patch_sf <- as.polygons(r_clean, dissolve = TRUE)
    patch_sf <- st_as_sf(patch_sf)
    final_freq <- as.data.frame(freq(r_clean, usena = FALSE))
    final_freq <- final_freq[,-1]
    colnames(final_freq) <- c("pid", "ncell")
    
    patch_sf <- patch_sf %>%
      left_join(final_freq, by = "pid") %>%
      mutate(patch_id = row_number())
    r_final <- rasterize(vect(patch_sf), r_clean, field = "patch_id")

    writeRaster(r_final, 
                file.path(patch_elim_rast_dir, paste0(fname, '.tif')),
                overwrite = TRUE)
    st_write(patch_sf,
             file.path(patch_elim_sf_dir, paste0(fname, '.shp')),
             delete_dsn = TRUE, quiet = TRUE)
  }
}

process_species_patches('Megaptera', year = 2003)
process_species_patches('Gadus', year = 2003)
process_species_patches('Hippoglossus', year = 2003)
