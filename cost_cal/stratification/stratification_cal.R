library(terra)
library(gsw)
library(stringr)

depth_lastringrdepth_layers <- list(
  c(0,5),c(5,10),c(10,25),c(25,50),c(50,100),c(100,150),c(150,200)
)
for(i in seq_along(depth_layers)){
  range_i <- depth_layers[[i]]
  dmin <- range_i[1]
  dmax <- range_i[2]
  temp <- rast(sprintf('T_layers/t_layer_%d_%dm_2003.tif',dmin,dmax))
  salt <- rast(sprintf('sal_layers/sal_layer_%d_%dm_2003.tif',dmin,dmax))
  z <- -(dmin + dmax)/2
  lonlat <- crds(temp,df = TRUE)
  lon <- lonlat[,1]
  lat <- lonlat[,2]
  p <- gsw_p_from_z(z,lat)
  SP <- values(salt)
  SA <- gsw_SA_from_SP(
    SP = SP,
    p = p,
    lon = lon,
    lat = lat
  )
  t <-values(temp)
  CT <- gsw_CT_from_t(
    SA = SA,
    t = t,
    p = p
  )
  rho <- gsw_rho(
    SA = SA,
    CT = CT,
    p = p
  )
  rho_rast <- temp
  values(rho_rast) <- rho
  out_file <- sprintf('rho_layers/rho_%d_%dm_2003.tif',dmin,dmax)
  writeRaster(rho_rast,out_file,overwrite = TRUE)
}
rho_files <- list.files(
  'rho_layers',
  pattern = '_2003\\.tif$',
  full.names = TRUE
)
depth_start <- as.numeric(str_extract(
  str_extract(rho_files, "rho_\\d+"), "\\d+"
))
rho_files_sorted <- rho_files[order(depth_start)]
rho_stack <- rast(rho_files_sorted)

depths <- c(2.5,7.5,17.5,37.5,75,125,175)
dz <- diff(depths)

drho_dz_list <- list()
for(k in 1:(nlyr(rho_stack)-1)){
  rho_upper <- rho_stack[[k]]
  rho_lower <- rho_stack[[k+1]]
  drho_dz <- (rho_lower - rho_upper)/-dz[k]
  drho_dz_list[[k]] <- drho_dz
}

g <- 9.81
rho0 <- 1025
N2_list <- list()
for(k in seq_along(drho_dz_list)){
  N2 <- -(g/rho0) * drho_dz_list[[k]]
  N2_list[[k]] <- N2
}

depth_names <- sapply(depth_layers, function(x) {
  paste0(x[1], "_", x[2], "m")
})

for (k in 1:7) {
  if (k == 1) {
    N2_layer_k <- N2_list[[1]]
  } else if (k == 7) {
    N2_layer_k <- N2_list[[6]]
  } else {
    N2_layer_k <- (N2_list[[k - 1]] + N2_list[[k]]) / 2
  }
  fname <- paste0("N2_", depth_names[k])

  tif_path <- file.path('n2_layers', paste0(fname, "_2003.tif"))
  writeRaster(N2_layer_k, tif_path, overwrite = TRUE)

  df <- as.data.frame(N2_layer_k, xy = TRUE, na.rm = TRUE)
  colnames(df) <- c("lon", "lat", "N2")
  csv_path <- file.path('n2_layers', paste0(fname, "_2003.csv"))
  write.csv(df, csv_path, row.names = FALSE)
}

