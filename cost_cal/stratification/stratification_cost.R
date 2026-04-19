library(dplyr)
library(readr)
library(ggplot2)
library(purrr)
library(stringr)
library(FNN)
library(data.table)
library(terra)
setwd('D:/Rwork/connectivity/n2_layers')

get_log_range <- function(n2_path) {
  df <- read.csv(n2_path)
  x_log <- log10(df$N2 * 1e4 + 1)
  return(c(min(x_log, na.rm = TRUE), max(x_log, na.rm = TRUE)))
}
n2_files <- list.files(pattern = '^N2_.*\\.csv$', full.names = TRUE) %>%
  sub("^\\./", "", .)
ranges <- map(n2_files, get_log_range)
global_xmin <- min(sapply(ranges, `[`, 1))
global_xmax <- max(sapply(ranges, `[`, 2))

scale_cost_global <- function(x, g_min, g_max){
  x_log <- log10(x * 1e4 + 1)
  scaled <- 1 + 99 * (x_log - g_min) / (g_max - g_min)
  return(scaled)
}

cost_stratification_cal <- function(n2_path,g_min, g_max){
  n2 <- read.csv(n2_path)
  fname <- basename(n2_path)
  depth_str <- str_extract(fname, '\\d+_\\d+m')
  n2_costs <- n2 %>%
    mutate(
      depth_layer=depth_str,
      costs=scale_cost_global(N2,g_min, g_max)
      )

  out_file<-paste0('cost_',n2_path)
  write.csv(n2_costs,out_file,row.names = FALSE)
}

walk(n2_files, ~cost_stratification_cal(.x, global_xmin, global_xmax))

cost_files <- list.files(pattern = '^cost_.*\\.csv$', full.names = TRUE)

for(f in cost_files){
  cc <- read.csv(f)
  fname <- basename(f)
  outname <- sub("\\.csv$", ".png", fname)
  p<-ggplot(cc, aes(x = lon, y = lat, fill = costs)) +
    geom_tile() +
    scale_fill_viridis_c(limits = c(0, 100)) +
    coord_fixed() +
    theme_minimal()
  ggsave(outname,p,width=8,height=6,dpi=600)
}
n2_tif_files <- list.files(pattern = '\\.tif$',full.names = TRUE)
for(f in n2_tif_files){
  n2 <- rast(f)
  fname <- basename(f)
  outname <- sub("\\.tif$", ".png", fname)
  png(outname, width = 800,height = 600,res = 300) 
  plot(n2,range = c(0,0.0006))
  dev.off()
}
