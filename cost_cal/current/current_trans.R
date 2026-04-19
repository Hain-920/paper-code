library(terra)
library(ggplot2)
input_dir<-'current_data/'
out_file<-'current_data_depth'
nc_files<-list.files(input_dir,pattern='\\.nc$',full.names = TRUE)
for(file in nc_files){
  fname<-basename(file)
  month_str <- sub(".*_(\\d{2})\\.nc$", "\\1", fname)
  current_u<-rast('EVEL_2003_01.nc')
  current_u1<-current_u[[1:18]]
  crs(current_u1)<-"EPSG:4326"
  ext(current_u1) <- c(-180, 179.9, -80, 90)
  current_u1<-flip(current_u1,'vertical')
  layer_list <- list(
    "0_5"     = current_u1[[1]],
    "5_10"    = mean(current_u1[[1:2]], na.rm = TRUE),
    "10_25"   = mean(current_u1[[2:3]], na.rm = TRUE),
    "25_50"   = mean(current_u1[[3:6]], na.rm = TRUE),
    "50_100"  = mean(current_u1[[6:11]], na.rm = TRUE),
    "100_150" = mean(current_u1[[11:15]], na.rm = TRUE),
    "150_200" = mean(current_u1[[15:18]], na.rm = TRUE)
  )
  for (depth_name in names(layer_list)) {
    out_file <- file.path(outdir,
                          paste0("current_u_", month_str, "_", depth_name, ".tif"))
    writeRaster(layer_list[[depth_name]], out_file, overwrite = TRUE)
  }
  
}
indir <- "output_current_u"
outdir <- "monthly_mean_by_depth"
dir.create(outdir, showWarnings = FALSE)
c_files <- list.files(indir, pattern = "^current_u_\\d{2}_.+\\.tif$", full.names = TRUE)

depth_layers <- unique(sub(".*_(\\d+_\\d+)\\.tif$", "\\1", basename(c_files)))

for (depth in depth_layers) {
  depth_files <- list.files(
    indir,
    pattern = paste0("_", depth, "\\.tif$"),
    full.names = TRUE
  )

  rs <- rast(depth_files)
  mean_r <- mean(rs, na.rm = TRUE)
  outfile <- file.path(outdir, paste0("current_u_", depth, ".tif"))
  writeRaster(mean_r, outfile, overwrite = TRUE)

}




