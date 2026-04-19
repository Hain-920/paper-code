library(sf)
library(data.table)
library(dplyr)
library(terra)
library(zeallot)
library(ggplot2)
library(rnaturalearth)
library(scales)
library(RANN)
setwd("D:/Rwork/connectivity")

calculate_species_connectivity <- function(species_name, beta_val = 0.01) {
  depth_layers <- c(
    "0_5","5_10","10_25","25_50",
    "50_100","100_150","150_200"
  )
  
  out_dir_sf <- file.path("patches", species_name, "trans_cost", "sf")
  out_dir_csv <- file.path("patches", species_name, "trans_cost", "csv")
  if(!dir.exists(out_dir_sf)) dir.create(out_dir_sf, recursive = TRUE)
  if(!dir.exists(out_dir_csv)) dir.create(out_dir_csv, recursive = TRUE)
  
  all_patch_results <- list()
  
  get_neighbor_layers <- function(layer){
    idx <- match(layer, depth_layers)
    up_layer <- NULL
    down_layer <- NULL
    
    if(idx > 1)
      up_layer <- depth_layers[idx-1]
    if(idx < length(depth_layers))
      down_layer <- depth_layers[idx+1]
    
    return(list(up = up_layer, down = down_layer))
  }
  
  # 投影计算
  project_csv <- function(csv_path, crs_target = 6933){
    dt <- fread(csv_path)
    sf_pts <- st_as_sf(dt, coords = c("lon", "lat"), crs = 4326)
    sf_proj <- st_transform(sf_pts, crs_target)
    coords <- st_coordinates(sf_proj)
    dt[, `:=`(x = coords[,1], y = coords[,2])]
    return(dt)
  }
  
  assign_patch_id <- function(cur_dt,str_dt, patch_sf){
    do_assign <- function(dt) {
      sf_obj <- st_as_sf(dt, coords = c("lon", "lat"), crs = 4326) %>% 
        st_transform(6933)
      
      nearest_idx <- st_nearest_feature(sf_obj, patch_sf_trans)
      
      dt$pid <- patch_sf_trans$pid[nearest_idx]
      return(dt)
    }
    cur_dt <- do_assign(cur_dt)
    str_dt <- do_assign(str_dt)
    return(list(cur = cur_dt, strat = str_dt))
  }
  
  project_dt <- function(data, crs_target = 6933){
    dt <- as.data.table(data)
    sf_pts <- st_as_sf(dt, coords = c("lon", "lat"), crs = 4326)
    sf_proj <- st_transform(sf_pts, crs_target)
    coords <- st_coordinates(sf_proj)
    dt[, `:=`(x = coords[,1], y = coords[,2])]
    return(dt)
  }

  # 计算单一格点转移成本
  calc_cell_cost <- function(cur_from, cur_to,
                             strat_from, strat_to,
                             gamma = 0.6){
    get_mid_depth <- function(depth_str){
      depth_str <- gsub("m", "", depth_str)
      nums <- as.numeric(strsplit(depth_str, "_")[[1]])
      mean(nums)
    }
    
    depth_from <- get_mid_depth(cur_from$depth_layer)
    depth_to   <- get_mid_depth(cur_to$depth_layer)
    
    dx <- cur_to$x - cur_from$x
    dy <- cur_to$y - cur_from$y
    dz <- depth_to - depth_from
    dist_m <- sqrt(dx^2 + dy^2 + dz^2)

    if (abs(dx) < 1e-3 && abs(dy) < 1e-3) {
      C_strat <- (strat_from + strat_to) / 2
      C_current <- (cur_from$costs + cur_to$costs) / 2
      C_total <- (C_strat + C_current) / 2
      
      return(C_total)
    } else {
      C_strat <- ((strat_from + strat_to)/2) * dist_m
      C_current <- ((cur_from$costs + cur_to$costs)/2) * dist_m
      
      mv_angle <- atan2(dy, dx)
      cur_angle <- atan2(cur_from$v, cur_from$u)
      cos_theta <- cos(mv_angle - cur_angle)
      f_theta <- 1 - gamma * cos_theta
      
      C_total <- (C_strat + f_theta * C_current) / 2
      
      return(C_total)
    }
  }
  
  # 找到与目标格点最近的格点
  nearest_cost <- function(query_xy, strat_dt){
    idx <- nn2(
      data = strat_dt[,c("x","y")],
      query = query_xy,
      k = 1
    )$nn.idx
    
    strat_dt$costs[idx]
  }
  
  calc_vertical_cost <- function(target_pid,
                                 cur_k, str_k,
                                 cur_down=NULL, cur_up=NULL,
                                 str_down=NULL, str_up=NULL,
                                 N0 = 50,
                                 beta = 0.01){
    get_mid_depth <- function(depth_str){
      depth_str <- gsub("m", "", depth_str)
      nums <- as.numeric(strsplit(depth_str, "_")[[1]])
      mean(nums)
    }
    pts <- cur_k[pid == target_pid]  
    if(nrow(pts) == 0) return(NULL)
    
    if(nrow(pts) > N0){
      pts <- pts[sample(.N, N0)]
    }   # 选取斑块内格点current成本
    
    cell_results <- lapply(1:nrow(pts), function(i){
      cell_i <- pts[i]
      q_xy <- matrix(c(cell_i$x, cell_i$y), 1, 2)
      strat_i <- nearest_cost(q_xy, str_k)   # 找到距离current格点最近strat数据,也就是格点strat成本
      depth_i <- get_mid_depth(cell_i$depth_layer)
      
      w_up <- NA
      w_down <- NA
      
      if(!is.null(cur_up)){
        idx_up <- nn2(data = cur_up[, c("x","y")], query = q_xy, k = 1)$nn.idx[1]
        cell_up <- cur_up[idx_up]
        strat_up <- nearest_cost(matrix(c(cell_up$x, cell_up$y), 1, 2), str_up)   # 上层strat成本 
        
        cost_up <- calc_cell_cost(cell_i, cell_up, strat_i, strat_up) 
        
        depth_up <- get_mid_depth(cell_up$depth_layer)
        dz_up <- abs(depth_up - depth_i)
        w_up <- cost_up * (1 + beta * dz_up)
      }
      
      if(!is.null(cur_down)){
        idx_down <- nn2(data = cur_down[, c("x","y")], query = q_xy, k = 1)$nn.idx[1]
        cell_down <- cur_down[idx_down]
        strat_down <- nearest_cost(matrix(c(cell_down$x, cell_down$y), 1, 2), str_down)
        
        cost_down <- calc_cell_cost(cell_i, cell_down, strat_i, strat_down)
        
        depth_down <- get_mid_depth(cell_down$depth_layer)
        dz_down <- abs(depth_down - depth_i)
        w_down <- cost_down * (1 + beta * dz_down)
      }
      
      if(is.na(w_up) & is.na(w_down)){
        return(c(connectivity=NA))
      }
      
      if (is.na(w_up)) {
        connectivity <- 1 / w_down 
      } else if (is.na(w_down)) {
        connectivity <- 1 / w_up
      } else {
        connectivity <- 1 / w_up + 1 / w_down 
      }
      return(c(connectivity = connectivity))
    })
    
    cell_df <- as.data.frame(do.call(rbind, cell_results))
    cell_df$connectivity <- as.numeric(cell_df$connectivity)
    
    return(list(
      pid = target_pid,
      mean_connectivity = mean(cell_df$connectivity, na.rm=TRUE)
    ))
  }
  
  calc_directional_cost <- function(target_pid,
                                    cur_k, str_k,
                                    cur_down=NULL, cur_up=NULL,
                                    str_down=NULL, str_up=NULL,
                                    N0 = 50){
    
    pts <- cur_k[pid == target_pid]  
    if(nrow(pts) == 0) return(NULL)
    
    if(nrow(pts) > N0){
      pts <- pts[sample(.N, N0)]
    }   
    
    cur_xy <- as.matrix(cur_k[, c("x", "y")])
    if(!is.null(cur_up)) up_xy <- as.matrix(cur_up[, c("x", "y")])
    if(!is.null(cur_down)) down_xy <- as.matrix(cur_down[, c("x", "y")])
    
    cell_results <- lapply(1:nrow(pts), function(i){
      cell_i <- pts[i]   # 起始格点current成本
      q_xy <- matrix(c(cell_i$x, cell_i$y), 1, 2)
      strat_i <- nearest_cost(q_xy, str_k)   # 起始strat成本
      
      candidates_cost <- numeric()
      candidates_dir <- character()   # 储存候选格点的数据
      
      idx_cur <- nn2(data = cur_xy, query = q_xy, k = 9)$nn.idx[1, 2:9]   # 起始格点最近8格点
      
      for(j in idx_cur){
        cell_j <- cur_k[j]
        strat_j <- nearest_cost(matrix(c(cell_j$x, cell_j$y), 1, 2), str_k)   # 最近格点成本值
        cst <- calc_cell_cost(cell_i, cell_j, strat_i, strat_j)
        
        if(!is.na(cst)){
          candidates_cost <- c(candidates_cost, cst)
          candidates_dir <- c(candidates_dir, "horizon")
        }
      }
      
      if(!is.null(cur_up)){
        idx_up <- nn2(data = up_xy, query = q_xy, k = 9)$nn.idx[1, 2:9]
        for(j in idx_up){
          cell_j <- cur_up[j]
          strat_j <- nearest_cost(matrix(c(cell_j$x, cell_j$y), 1, 2), str_up)
          cst <- calc_cell_cost(cell_i, cell_j, strat_i, strat_j)
          
          if(!is.na(cst)){
            candidates_cost <- c(candidates_cost, cst)
            candidates_dir <- c(candidates_dir, "up")
          }
        }
      }
      
      if(!is.null(cur_down)){
        idx_down <- nn2(data = down_xy, query = q_xy, k = 9)$nn.idx[1, 2:9]
        for(j in idx_down){
          cell_j <- cur_down[j]
          strat_j <- nearest_cost(matrix(c(cell_j$x, cell_j$y), 1, 2), str_down)
          cst <- calc_cell_cost(cell_i, cell_j, strat_i, strat_j)
          
          if(!is.na(cst)){
            candidates_cost <- c(candidates_cost, cst)
            candidates_dir <- c(candidates_dir, "down")
          }
        }
      }
      
      if(length(candidates_cost) == 0){
        return(c(cost=NA, dir=NA))
      }
      
      min_idx <- which.min(candidates_cost)
      min_cost <- candidates_cost[min_idx]
      best_dir <- candidates_dir[min_idx]
      
      return(c(cost=min_cost, dir=best_dir))
    })
    
    cell_df <- as.data.frame(do.call(rbind, cell_results))
    cell_df$cost <- as.numeric(cell_df$cost)
    
    return(list(
      pid               = target_pid,          
      mean_min_cost     = mean(cell_df$cost, na.rm=TRUE),
      current_layer_rat = mean(cell_df$dir == "horizon", na.rm=TRUE),
      up_layer_rat      = mean(cell_df$dir == "up", na.rm=TRUE),
      down_layer_rat    = mean(cell_df$dir == "down", na.rm=TRUE)
    ))
  }
  
  for(layer in depth_layers){
    cat("Processing layer:", layer, "\n")
    neigh <- get_neighbor_layers(layer)
    
    up_layer <- neigh$up
    down_layer <- neigh$down
    cat("Upper layer:", up_layer, "\n")
    cat("Lower layer:", down_layer, "\n")
    
    patch_shp_path <- file.path("patches", species_name, "patches_elim", "sf", paste0("patch_", layer, "m_2003.shp"))
    if(!file.exists(patch_shp_path)) {
      next
    }
    
    # 当前层数据
    cur_dt <- fread(paste0("merged_current/cost_current_",layer,".csv"))
    str_dt <- fread(paste0("n2_layers/cost_N2_",layer,"m_2003.csv"))
    cur_k <- project_dt(cur_dt)
    str_k <- project_dt(str_dt)
    
    # 斑块数据
    patch_sf <- st_read(patch_shp_path, quiet=TRUE)
    patch_sf_trans <- st_transform(patch_sf, 6933)
    c(cur_k, str_k) %<-% assign_patch_id(cur_k, str_k, patch_sf_trans)

    # 上层
    cur_up <- NULL; str_up <- NULL
    if(!is.null(up_layer)){
      cat("Loading upper layer:", up_layer, "\n")
      cur_up <- project_csv(
        paste0("merged_current/cost_current_",up_layer,".csv")
      )
      str_up <- project_csv(
        paste0("n2_layers/cost_N2_",up_layer,"m_2003.csv")
      )
    } 
    
    # 下层
    cur_down <- NULL; str_down <- NULL
    if(!is.null(down_layer)){
      cat("Loading lower layer:", down_layer, "\n")
      cur_down <- project_csv(
        paste0("merged_current/cost_current_",down_layer,".csv")
      )
      str_down <- project_csv(
        paste0("n2_layers/cost_N2_",down_layer,"m_2003.csv")
      )
    } 
    
    pids <- unique(cur_k$pid)
    cat("Calculating vertical cost\n")
    vertical_results <- lapply(pids, function(p){
      set.seed(123)
      calc_vertical_cost(
        target_pid = p,
        cur_k      = cur_k, 
        str_k      = str_k,
        cur_down   = cur_down, 
        cur_up     = cur_up,
        str_down   = str_down, 
        str_up     = str_up,
        beta       = beta_val
      )
    })
    vertical_df <- do.call(rbind,
                           lapply(vertical_results, as.data.frame))
    
    cat("Calculating directional cost\n")
    directional_results <- lapply(pids, function(p){
      calc_directional_cost(
        target_pid = p,
        cur_k      = cur_k,
        str_k      = str_k,
        cur_down   = cur_down,
        cur_up     = cur_up,
        str_down   = str_down,
        str_up     = str_up
      )
    })
    directional_df <- do.call(rbind,
                              lapply(directional_results, as.data.frame))
    
    res_df <- vertical_df %>%
      left_join(directional_df, by = "pid")
    
    patch_dt <- patch_sf %>%
      left_join(res_df, by = "pid") %>%
      mutate(
        connectivity_log = log10(mean_connectivity),
        dominant_dir = case_when(
          current_layer_rat >= up_layer_rat &
            current_layer_rat >= down_layer_rat ~ "horizon",
          up_layer_rat > current_layer_rat &
            up_layer_rat >= down_layer_rat ~ "up",
          down_layer_rat > current_layer_rat &
            down_layer_rat > up_layer_rat ~ "down",
          TRUE ~ "unknown"),
        signed_cost_directional = case_when(
          dominant_dir == "up"      ~ mean_min_cost,
          dominant_dir == "down"    ~ -mean_min_cost,
          dominant_dir == "horizon" ~ 0,
          TRUE ~ NA_real_),
        depth_layer = layer
      )
    out_shp_path <- file.path(out_dir_sf, paste0("patch_cost_", layer, "m_2003.gpkg"))
    st_write(patch_dt, 
             out_shp_path, 
             append = FALSE, quiet = TRUE)
    
    all_patch_results[[layer]] <- patch_dt
  }
  patch_dt_all <- bind_rows(all_patch_results)
  patch_dt_all$depth_layer <- factor(
    patch_dt_all$depth_layer, 
    levels = c("0_5", "5_10", "10_25", "25_50", "50_100", "100_150", "150_200")
  )
  out_all_shp_path <- file.path(out_dir_sf, 'all_layer_cost.gpkg')
  st_write(patch_dt_all, 
           out_all_shp_path, 
           append = FALSE, quiet = TRUE)
  return(patch_dt_all)
}
h_patch_dt <- calculate_species_connectivity('Hippoglossus',beta_val = 0.01)
g_patch_dt <- calculate_species_connectivity('Gadus',beta_val = 0.01)
m_patch_dt <- calculate_species_connectivity('Megaptera',beta_val = 0.01)

