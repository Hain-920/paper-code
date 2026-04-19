library(ggtern)
library(ggplot2)
library(dplyr)

world <- ne_countries(scale = "medium", returnclass = "sf")
plot_dominant_direction <- function(patch_data, species_name,world_data = world) {
  p_facet_dominant <- ggplot() +
    geom_sf(data = world_data, fill = "gray85", color = "white", size = 0.2) +
    geom_sf(data = patch_data, aes(fill = dominant_dir), color = 'black', size = 0.1) + 
    scale_fill_manual(
      breaks = c("up", "down", "horizon", "unknown"),
      values = c("up"      = "#d73027",  
                 "down"    = "#4575b4",  
                 "horizon" = "#ffffbf",  
                 "unknown" = "grey50"),
      labels = c("Up", "Down", "Horizontal", "Unknown"),
      name = "Dominant Direction"
    ) +
    facet_wrap(~ depth_layer, ncol = 3) + 
    theme_minimal() +
    coord_sf(xlim = st_bbox(patch_data)[c(1,3)], ylim = st_bbox(patch_data)[c(2,4)]) + 
    labs(title = bquote(italic(.(species_name)) ~ "Dominant Migration Direction")) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      panel.background = element_rect(fill = "#e0f3f8", color = NA),
      panel.grid.major = element_blank(),
      strip.background = element_rect(fill = "gray30"),
      strip.text = element_text(color = "white", face = "bold"),
      legend.position = "bottom"
    ) 
  out_name <- paste0(gsub(" ", "_", species_name), "_dominant_direction.png")
  ggsave(filename = out_name, plot = p_facet_dominant, width = 8, height = 8, dpi = 600)
  return(p_facet_dominant)
}
plot_dominant_direction(g_patch_dt,'Gadus',world)
plot_dominant_direction(h_patch_dt,'Hippoglossus',world)
plot_dominant_direction(m_patch_dt,'Megaptera',world)

plot_vertical_connectivity <- function(patch_data, species_name, world_data = world) {
  p_facet_cost <- ggplot() +
    geom_sf(data = world_data, fill = "gray85", color = "white", size = 0.2) +
    geom_sf(data = patch_data,
            aes(fill = connectivity_log),
            color = 'grey40', linewidth = 0.2, alpha = 0.9) + 
    scale_fill_viridis_c(
      option = "turbo", 
      name = "Vertical\nConnectivity (log)", 
      na.value = "transparent"
    ) +
    facet_wrap(~ depth_layer, ncol = 3) + 
    coord_sf(
      xlim = st_bbox(patch_data)[c(1,3)],
      ylim = st_bbox(patch_data)[c(2,4)]
    ) + 
    theme_bw() +
    labs(title = bquote(italic(.(species_name)))) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      panel.background = element_rect(fill = "#e0f3f8", color = NA),
      panel.grid.major = element_line(color = "grey90", linetype = "dashed"),
      strip.background = element_rect(fill = "gray30"),
      strip.text = element_text(color = "white", face = "bold"),
      legend.position = "bottom",
      legend.key.width = unit(3, "cm"),     
      legend.key.height = unit(0.5, "cm"),  
      legend.title = element_text(vjust = 0.8), 
      legend.text = element_text(size = 8)
    )

  out_filename <- paste0(gsub(" ", "_", species_name), '_vertical_connectivity.png')
  ggsave(filename = out_filename, plot = p_facet_cost, height = 8, width = 8, dpi = 600)
  
  return(p_facet_cost)
}
plot_vertical_connectivity(h_patch_dt,'Hippoglossus',world)
plot_vertical_connectivity(g_patch_dt,'Gadus',world)
plot_vertical_connectivity(m_patch_dt,'Megaptera',world)

plot_habitat_fragmentation <- function(raw_patch_df, species_name) {
  stats_df <- raw_patch_df %>%
    group_by(depth_layer) %>%
    summarise(
      cell_count = sum(ncell, na.rm = TRUE),
      patch_count = n_distinct(pid) 
    ) %>%
    ungroup()
  
  stats_df$depth_layer <- factor(stats_df$depth_layer, 
                                 levels = c("0_5", "5_10", "10_25", "25_50", "50_100", "100_150", "150_200"))
  stats_df <- stats_df[!is.na(stats_df$depth_layer), ]
  scale_factor <- max(stats_df$patch_count, na.rm = TRUE) / max(stats_df$cell_count, na.rm = TRUE)
  
  p_frag <- ggplot(stats_df, aes(x = depth_layer)) +
    geom_col(aes(y = patch_count), fill = "#74add1", alpha = 0.8, width = 0.6) +

    geom_line(aes(y = cell_count * scale_factor, group = 1), color = "#d73027", linewidth = 1.2) +
    geom_point(aes(y = cell_count * scale_factor), color = "#d73027", size = 3, shape = 21, fill = "white", stroke = 1.5) +

    scale_y_continuous(
      name = "Number of Patches (Fragmentation Degree)",
      sec.axis = sec_axis(~ . / scale_factor, name = "Total Cell Count (Habitat Capacity)")
    ) +
    
    theme_bw() +
    labs(
      title = bquote(italic(.(species_name)) ~ "Habitat Fragmentation across Depths"),
      x = "Depth Layer (m)"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title.y.left = element_text(color = "#313695", face = "bold", margin = ggplot2::margin(r = 10)),
      axis.text.y.left = element_text(color = "#313695", face = "bold"),
      axis.title.y.right = element_text(color = "#d73027", face = "bold", margin = ggplot2::margin(l = 10)),
      axis.text.y.right = element_text(color = "#d73027", face = "bold"),
      
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      panel.grid.minor = element_blank()
    )
  out_name <- paste0(gsub(" ", "_", species_name), "_fragmentation_swapped.png")
  ggsave(out_name, p_frag, width = 8, height = 5, dpi = 600)
  
  return(p_frag)
}
plot_habitat_fragmentation(h_patch_dt,'Hippoglossus')
plot_habitat_fragmentation(m_patch_dt,'Megaptera')
plot_habitat_fragmentation(g_patch_dt,'Gadus')


