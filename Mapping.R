pacman::p_load(tidyverse,linelist,rio,here,sp,sf,raster,tmap,foreign,ncdf4,data.table,terra,rgdal,MODISTools,tmap,ggplot2,viridis,ggnewscale,ggspatial,ggtext,cowplot,ggpubr)

# Import the watershed boundary shapefile
basin_gsl <- read_sf('~/basin_gsl/basin_gsl.shp')
alfalfa_wta <- import('~/alfalfa_wta.csv') # Note: alfalfa_wta dataframe can also be obtained using code in Water WTA.R

####### Plot Figs 1 and 2 #######
alfalfa_wta_fallow <- alfalfa_wta %>% rename(value = wta_Fallow) %>% mutate(color = NA) %>% dplyr::select(X,Y,value,color)
alfalfa_wta_grain <- alfalfa_wta %>% rename(value = wta_Grain) %>% mutate(color = NA) %>% dplyr::select(X,Y,value,color)
alfalfa_wta_hay <- alfalfa_wta %>% rename(value = wta_Hay) %>% mutate(color = NA) %>% dplyr::select(X,Y,value,color)

wta_limits <- range(c(alfalfa_wta_fallow$value,alfalfa_wta_grain$value,alfalfa_wta_hay$value),na.rm = TRUE)
wta_breaks <- seq(2500,6500,by=500)
wta_colors <- c('#01665e','#35978f','#80cdc1','#c7eae5','#f6e8c3','#dfc27d','#bf812d','#8c510a')  

alfalfa_cost_fallow <- alfalfa_wta %>% rename(value = cost_Fallow) %>% mutate(color = NA) %>% dplyr::select(X,Y,value,color)
alfalfa_cost_grain <- alfalfa_wta %>% rename(value = cost_Grain) %>% mutate(color = NA) %>% dplyr::select(X,Y,value,color)
alfalfa_cost_hay <- alfalfa_wta %>% rename(value = cost_Hay) %>% mutate(color = NA) %>% dplyr::select(X,Y,value,color)

cost_limits <- range(c(alfalfa_cost_fallow$value,alfalfa_cost_grain$value,alfalfa_cost_hay$value),na.rm = TRUE)
cost_breaks <- pretty(cost_limits)
cost_colors <- c('#01665e','#35978f','#80cdc1','#c7eae5','#f6e8c3','#dfc27d','#bf812d','#8c510a')  

# Define createMap function
createMap <- function(.data, .name, .color, .limits, .breaks, maptitle){
  ggplot() +
    geom_tile(data = .data, aes(x = X, y = Y, fill = value)) +
    scale_fill_gradientn(name = .name, colors = .color, limits = .limits, breaks = .breaks) +
    geom_sf(data = basin_gsl, fill = NA, color = "gray30") +
    geom_sf_text(data = basin_gsl, aes(label = basin_name), color = "gray30") + 
    annotation_scale(location = "bl", width_hint = 0.2, text_cex = 0.75, bar_cols = c("gray30")) +
    theme_minimal() +
    labs(title = maptitle, fill = "Value") +
    theme(
      plot.title = element_text(size = 12, hjust = 0.5, vjust = 1),  # Adjust the size, face, and alignment
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    ) +
    coord_sf()
}

map_wta_fallow <- createMap(.data = alfalfa_wta_fallow, maptitle = 'a. To Fallow', 
                            .name = 'WTA ($/ha)', .color = wta_colors, .limits = wta_limits, .breaks = wta_breaks)
map_wta_grain <- createMap(.data = alfalfa_wta_grain, maptitle = 'b. To Spring Grains', 
                           .name = 'WTA ($/ha)', .color = wta_colors, .limits = wta_limits, .breaks = wta_breaks)
map_wta_hay <- createMap(.data = alfalfa_wta_hay, maptitle = 'c. To Other Hays',
                         .name = 'WTA ($/ha)', .color = wta_colors, .limits = wta_limits, .breaks = wta_breaks)
map_wta_combined <- ggarrange(map_wta_fallow,map_wta_grain,map_wta_hay,nrow=1,ncol=3,common.legend = T,legend = "right") # Fig 2

map_cost_fallow <- createMap(.data = alfalfa_cost_fallow, maptitle = 'a. To Fallow',
                             .name = 'Cost ($/m3)', .color = cost_colors, .limits = cost_limits, .breaks = cost_breaks)
map_cost_grain <- createMap(.data = alfalfa_cost_grain, maptitle = 'b. To Spring Grains',
                            .name = 'Cost ($/m3)', .color = cost_colors, .limits = cost_limits, .breaks = cost_breaks)
map_cost_hay <- createMap(.data = alfalfa_cost_hay, maptitle = 'c. To Other Hays',
                          .name = 'Cost ($/m3)', .color = cost_colors, .limits = cost_limits, .breaks = cost_breaks)
map_cost_combined <- ggarrange(map_cost_fallow,map_cost_grain,map_cost_hay,nrow=1,ncol=3,common.legend = T, legend = "right") # Fig 2

####### Plot Figs 4 and 5 #######
alfalfa_wta_basin <- import('~\alfalfa_wta_basin.csv') 
alfalfa_wta_county <- import('~\alfalfa_wta_county.csv') # Note: alfalfa_wta_basin and alfalfa_wta_county can also obtained using code in Water WTA.R

# Define createMap_a function
createMap_a <- function(.value, .data, maptitle){
  ggplot() +
    geom_tile(data = .data, aes(x = X, y = Y, fill = factor({{.value}}))) +
    scale_fill_manual(name = 'Projected Status', values = c("1" = "#1a9850","0" = "#d73027"), labels = c("1" = "Enrolled","0" = "Not enrolled")) +
    geom_sf(data = basin_gsl, fill = NA, color = "gray30") +
    geom_sf_text(data = basin_gsl, aes(label = basin_name), color = "gray30") + 
    annotation_scale(location = "bl", width_hint = 0.2, text_cex = 0.75, bar_cols = c("gray30")) +
    theme_minimal() +
    labs(title = maptitle, fill = "Value") +
    theme(
      plot.title = element_markdown(size = 12, hjust = 0.5, vjust = 1),  # Adjust the size, face, and alignment
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    ) +
    coord_sf()
}

### Watershed-level payment system: Fig 4 ###
map_baseline <- createMap_a(maptitle = '', .value = enroll, .data = alfalfa_wta_basin)
map_baseline_lower <- createMap_a(maptitle = '', .value = enroll_upper, .data = alfalfa_wta_basin)
map_baseline_upper <- createMap_a(maptitle = '', .value = enroll_lower, .data = alfalfa_wta_basin)
map_conserved <- createMap_a(maptitle = '', .value = enroll_a, .data = alfalfa_wta_basin)
map_conserved_lower <- createMap_a(maptitle = '', .value = enroll_a_upper, .data = alfalfa_wta_basin)
map_conserved_upper <- createMap_a(maptitle = '', .value = enroll_a_lower, .data = alfalfa_wta_basin)

# Arrange the plots in a 2x3 grid
map_combined_basin <- ggarrange(map_baseline, map_baseline_lower,map_baseline_upper, map_conserved,map_conserved_lower,map_conserved_upper, 
                                ncol = 3, nrow = 2, common.legend = T, legend = "right", labels = "AUTO")

# Annotate the figure with row and column labels
annotated_plot_basin <- annotate_figure(map_combined_basin,
                                        top = text_grob(c(''), size = 12, face = "bold"),
                                        left = text_grob(c(''), size = 12, face = "bold", rot = 90))

ggdraw() +
  draw_plot(annotated_plot_basin, hjust = 0, vjust = 0, scale = 1) +
  draw_label("Baseline Scenario", x = 0.02, y = 0.75, angle = 90, size = 14, fontface = "bold") +
  draw_label("Conservative Scenario", x = 0.02, y = 0.25, angle = 90, size = 14, fontface = "bold") +
  draw_label("Projection", x = 0.14, y = 0.98, size = 14, fontface = "bold", vjust = 1) +
  draw_label("Lower Bound of 95% CI", x = 0.48, y = 0.98, size = 14, fontface = "bold", vjust = 1) +
  draw_label("Upper Bound of 95% CI", x = 0.78, y = 0.98, size = 14, fontface = "bold", vjust = 1)

### County-level payment system: Fig 5 ###
map_baseline <- createMap_a(maptitle = '', .value = enroll, .data = alfalfa_wta_county)
map_baseline_lower <- createMap_a(maptitle = '', .value = enroll_upper, .data = alfalfa_wta_county)
map_baseline_upper <- createMap_a(maptitle = '', .value = enroll_lower, .data = alfalfa_wta_county)
map_conserved <- createMap_a(maptitle = '', .value = enroll_a, .data = alfalfa_wta_county)
map_conserved_lower <- createMap_a(maptitle = '', .value = enroll_a_upper, .data = alfalfa_wta_county)
map_conserved_upper <- createMap_a(maptitle = '', .value = enroll_a_lower, .data = alfalfa_wta_county)

# Arrange the plots in a 2x3 grid
map_combined_county <- ggarrange(map_baseline, map_baseline_lower,map_baseline_upper, map_conserved,map_conserved_lower,map_conserved_upper, 
                                 ncol = 3, nrow = 2, common.legend = T, legend = "right", labels = "AUTO")

# Annotate the figure with row and column labels
annotated_plot_county <- annotate_figure(map_combined_county,
                                         top = text_grob(c(''), size = 12, face = "bold"),
                                         left = text_grob(c(''), size = 12, face = "bold", rot = 90))

ggdraw() +
  draw_plot(annotated_plot_county, hjust = 0, vjust = 0, scale = 1) +
  draw_label("Baseline Scenario", x = 0.02, y = 0.75, angle = 90, size = 14, fontface = "bold") +
  draw_label("Conservative Scenario", x = 0.02, y = 0.25, angle = 90, size = 14, fontface = "bold") +
  draw_label("Projection", x = 0.14, y = 0.98, size = 14, fontface = "bold", vjust = 1) +
  draw_label("Lower Bound of 95% CI", x = 0.48, y = 0.98, size = 14, fontface = "bold", vjust = 1) +
  draw_label("Upper Bound of 95% CI", x = 0.78, y = 0.98, size = 14, fontface = "bold", vjust = 1)

