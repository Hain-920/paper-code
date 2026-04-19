library(geosphere)
library(dplyr)
library(ggplot2)
library(purrr)
library(stringr)
library(FNN)
library(data.table)
library(terra)
setwd('D:/Rwork/connectivity/merged_current')

# 将时间成本范围映射到1-100
get_days_log_range <- function(current_path){
  df <- read.csv(current_path)[,-1]
  dist_ns_km <- 111.32
  dist_ew_km <- pmax(111.32*cos(df$lat*pi/180),5)
  current_angle <- atan2(df$v, df$u)
  
  effective_dist_km <- ifelse(
    df$speed < 0.01,
    pmax(dist_ew_km, dist_ns_km),
    ifelse(
      abs(current_angle) < pi/4 | abs(current_angle) > 3*pi/4,
      dist_ew_km/abs(cos(current_angle)),
      dist_ns_km/abs(sin(current_angle))
    )
  )
  cost_precise_days <- effective_dist_km/(pmax(df$speed,0.05)*86.4)
  x_log <- log10(cost_precise_days + 1)
  
  return(c(min(x_log,na.rm=TRUE),max(x_log,na.rm=TRUE)))
}
current_files <- list.files(pattern = '^current_.*\\.csv$', full.names = TRUE) %>%
  sub("^\\./", "", .)
ranges <- purrr::map(current_files,get_days_log_range)

global_xmin <- min(sapply(ranges,`[`,1))
global_xmax <- max(sapply(ranges,`[`,2))

scale_cost_days <- function(x,g_min,g_max){
  x_log <- log10(x+1)
  scaled <- 1 + 99*(x_log - g_min)/(g_max - g_min)
  
  return(scaled)
}
cost_current_cal<-function(current_path){
  current<-read.csv(current_path)[,-1]
  fname<-basename(current_path)
  depth_str<-str_extract(fname,'\\d+_\\d+')
  depth_str<-depth_str%>%paste0('m')
  
  current<-current%>%
    mutate(
      speed = ifelse(speed <= 0, NA, speed),
      safe_speed = ifelse(is.na(speed), NA, pmax(speed,0.05)),
      depth_layer=depth_str,
      dist_ns_km=111.32,
      dist_ew_km=pmax(111.32*cos(lat*pi/180),5),
      current_angle=atan2(v, u),
      effective_dist_km = case_when(
        speed < 0.01 ~ pmax(dist_ew_km, dist_ns_km),
        abs(current_angle) < pi/4 | abs(current_angle) > 3*pi/4 ~ dist_ew_km / abs(cos(current_angle)),
        TRUE ~ dist_ns_km / abs(sin(current_angle))
      ),
      cost_precise_days = effective_dist_km/(safe_speed*86.4),
      costs = scale_cost_days(cost_precise_days,global_xmin,global_xmax)
    ) 
  
  out_file<-paste0('cost_',current_path)
  write.csv(current,out_file,row.names = FALSE)
}
walk(current_files,cost_current_cal)

cost_files <- list.files(pattern = '^cost_.*\\.csv$', full.names = TRUE)

for(f in cost_files){
  cc <- read.csv(f)
  fname <- basename(f)
  outname <- sub("\\.csv$", ".png", fname)
  p<-ggplot(cc, aes(x = lon, y = lat, fill = costs)) +
    geom_tile() +
    scale_fill_viridis_c(na.value = 'white') +
    coord_fixed() +
    theme_minimal()
  ggsave(outname,p,width=8,height=6,dpi=600)
}

cost_files<-sub("^\\./", "",list.files(pattern='cost.*\\.csv$', full.names = TRUE))
data_list<-lapply(cost_files,function(f){
  read.csv(f,stringsAsFactors = FALSE)
})
merged_costs<-do.call(rbind,data_list)
write.csv(merged_costs,'all_current_costs.csv',row.names = FALSE)

