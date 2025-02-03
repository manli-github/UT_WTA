pacman::p_load(tidyverse,linelist,rio,here,sp,sf,raster,tmap,foreign,ncdf4,data.table,terra,rgdal,MODISTools,tmap,ggplot2,viridis,ggnewscale,ggspatial,ggtext,cowplot,ggpubr)

# Import the watershed boundary shapefile
basin_gsl <- read_sf('~/basin_gsl/basin_gsl.shp')
lake <- read_sf('~/UtahLakesNHD/LakesNHDHighRes.shp')
study_area <- read_sf('~/study_area.shp')
ut <- read_sf('~/Utah/Utah.shp')
usa <- read_sf('~/gadm41_USA_shp/gadm41_USA_1.shp')
cdl2022 <- raster('~/cdl2022.tif')

alfalfa_wta <- import('~/alfalfa_wta.csv') # Note: alfalfa_wta dataframe can also be obtained using code in Water WTA.R

# Study area map: Fig. 1
lake_gsl <- lake %>% filter(GNIS_Name == 'Great Salt Lake') %>% 
  dplyr::select(GNIS_Name)
lake_gsl<- st_transform(lake_gsl, crs(study_area))
ut <- ut %>% filter(STATE=='Utah') %>% mutate(STATE='UTAH')
ut <- st_transform(ut, crs(study_area))

usa <- usa %>% filter(NAME_1 != 'Alaska' & NAME_1 != 'Hawaii')
usa <- st_union(usa)
usa <- st_transform(usa, crs(study_area))
county <- st_transform(county, crs(study_area))

crs(cdl2022) <- crs(study_area)
cdl_gsl <- crop(cdl2022,extent(lake_gsl))
cdl_gsl <- mask(cdl_gsl,lake_gsl)
all_ids <- unique(values(cdl_gsl))
reclass_vector <- rep(NA, length(all_ids))
reclass_vector[3] <- 111 # open water ID = 111
rclmat <- cbind(all_ids, reclass_vector)
cdl_water <- reclassify(cdl_gsl, rclmat)
cdl_water_df <- as.data.frame(cdl_water, xy = TRUE) %>% 
  rename (value = Class_Names) %>% 
  drop_na(value) %>% 
  mutate(color = ifelse(value == 111, "#3182bd", NA)) 

map_ut <- ggplot() +
  geom_sf(data = lake_gsl, fill = "#c6dbef", color = "#9ecae1") +
  geom_tile(data = cdl_water_df, aes(x = x, y = y, fill = color), na.rm = TRUE) +
  scale_fill_manual(name = 'Great Salt Lake', values = '#3182bd', labels = NULL, na.value = NA) +
  new_scale_fill() +  
  geom_sf(data = county, aes(color = "County Boundary"), fill = NA, size = 0.5, show.legend = TRUE) +
  geom_sf(data = basin_gsl, aes(color = "Watershed Boundary"), fill = scales::alpha('gray80', 0.5), size = 1.5, show.legend = TRUE) +
  geom_sf_text(data = county, aes(label = NAME), color = "gray30", size = 2.5) + 
  annotation_scale(location = "bl", width_hint = 0.2, text_cex = 0.75, bar_cols = c("gray30","white")) +
  scale_color_manual(name = "Boundaries", 
                     values = c("Watershed Boundary" = "red", "County Boundary" = "gray50"),
                     labels = c("County Boundary", "Watershed Boundary")) +
  annotation_north_arrow(
    location = "tl", which_north = "true",
    style = north_arrow_nautical(
      fill = c("grey40", "white"),
      line_col = "grey20")) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),       # Remove x-axis text
    axis.text.y = element_blank(),       # Remove y-axis text
    axis.ticks = element_blank(),        # Remove axis ticks
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    legend.title = element_text(size = 14),         
    legend.text = element_text(size = 14), 
    legend.key.size = unit(1.2, 'cm') 
  ) +
  coord_sf()

inset_map_a <- ggplot() +
  geom_sf(data = usa, fill = NA, color = "gray40") +
  geom_sf(data = ut, fill = "gray70", color = "gray70",alpha = 0.5) +
  geom_sf(data = basin_gsl, fill = NA, color = "red", size = 1.5, ) +
  theme_void() 

map_studyarea <- ggdraw() +
  draw_plot(map_ut, hjust = 0, vjust = 0, scale = 1) +
  draw_plot(inset_map_a, halign = 1, valign = 1, scale = .25)
ggsave('~/map_studyarea.jpeg',plot=map_studyarea,
       width=16,height=10,dpi=300)

####### Plot Figs 2 and 3 #######
alfalfa_wta_fallow <- alfalfa_wta %>% rename(value = wta_Fallow) %>% mutate(color = NA) %>% dplyr::select(X,Y,value,color)
alfalfa_wta_grain <- alfalfa_wta %>% rename(value = wta_Grain) %>% mutate(color = NA) %>% dplyr::select(X,Y,value,color)
alfalfa_wta_hay <- alfalfa_wta %>% rename(value = wta_Hay) %>% mutate(color = NA) %>% dplyr::select(X,Y,value,color)

alfalfa_cost_fallow <- alfalfa_wta %>% rename(value = cost_Fallow) %>% mutate(color = NA) %>% dplyr::select(X,Y,value,color)
alfalfa_cost_grain <- alfalfa_wta %>% rename(value = cost_Grain) %>% mutate(color = NA) %>% dplyr::select(X,Y,value,color)
alfalfa_cost_hay <- alfalfa_wta %>% rename(value = cost_Hay) %>% mutate(color = NA) %>% dplyr::select(X,Y,value,color)

### WTA ###
# Define createMap function
createMap <- function(.data, .name, .color, .limits, .breaks, maptitle){
  ggplot() +
    geom_sf(data = lake_gsl, fill = "#c6dbef", color = "#9ecae1") +
    geom_tile(data = cdl_water_df, aes(x = x, y = y, fill = color), na.rm = TRUE, show.legend = FALSE) +
    scale_fill_manual(name = '', values = '#3182bd', labels = NULL, na.value = NA) +
    new_scale_fill() +  
    geom_tile(data = .data, aes(x = X, y = Y, fill = value)) +
    scale_fill_gradientn(name = .name, colors = .color, limits = .limits, breaks = .breaks) +
    geom_sf(data = basin_gsl, fill = NA, color = "gray30") +
    geom_sf_text(data = basin_gsl, aes(label = basin_name), color = "gray30") + 
    annotation_scale(location = "br", width_hint = 0.2, text_cex = 0.75, bar_cols = c("gray30","white")) +
    theme_minimal() +
    labs(title = maptitle, fill = "Value") +
    theme(
      plot.title = element_text(size = 16, hjust = 0.5, vjust = 1),  # Adjust the size, face, and alignment
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 14),
      legend.key.size = unit(1.2, 'cm')
    ) +
    coord_sf()
}

createMap_noscale <- function(.data, .name, .limits, .breaks, .value, maptitle){
  ggplot() +  
    geom_sf(data = lake_gsl, fill = "#c6dbef", color = "#9ecae1") +
    geom_tile(data = cdl_water_df, aes(x = x, y = y, fill = color), na.rm = TRUE, show.legend = FALSE) +
    scale_fill_manual(name = '', values = '#3182bd', labels = NULL, na.value = NA) +
    new_scale_fill() +  
    geom_tile(data = .data, aes(x = X, y = Y, fill = {{.value}})) +
    scale_fill_viridis_c(name = .name, limits = .limits, breaks = .breaks) +
    geom_sf(data = basin_gsl, fill = NA, color = "gray30") +
    geom_sf_text(data = basin_gsl, aes(label = basin_name), color = "gray30") + 
    theme_minimal() +
    labs(title = maptitle, fill = "Value") +
    theme(
      plot.title = element_text(size = 16, hjust = 0.5, vjust = 1,face = "bold"),  # Adjust the size, face, and alignment
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 14),
      legend.key.size = unit(1.2, 'cm')
    ) +
    coord_sf()
}

wta_limits_fallow <- range(c(alfalfa_wta_fallow$value),na.rm = TRUE)
wta_breaks_fallow <- pretty(wta_limits_fallow)
map_wta_fallow <- createMap_noscale(.data = alfalfa_wta_fallow, maptitle = '(a) To Fallow', 
                            .name = 'WTA ($/ha)', .limits = wta_limits_fallow, .breaks = wta_breaks_fallow, .value = value)
wta_limits_grain <- range(c(alfalfa_wta_grain$value),na.rm = TRUE)
wta_breaks_grain <- pretty(wta_limits_grain)
map_wta_grain <- createMap_noscale(.data = alfalfa_wta_grain, maptitle = '(b) To Spring Grains', 
                                   .name = 'WTA ($/ha)', .limits = wta_limits_grain, .breaks = wta_breaks_grain,  .value = value)
wta_limits_hay <- range(c(alfalfa_wta_hay$value),na.rm = TRUE)
wta_breaks_hay <- pretty(wta_limits_hay)
map_wta_hay <- createMap(.data = alfalfa_wta_hay, maptitle = '(c\u200B) To Other Hays', 
                         .name = 'WTA ($/ha)', .limits = wta_limits_hay, .breaks = wta_breaks_hay, .value = value)
map_wta_combined <- ggarrange(map_wta_fallow,map_wta_grain,map_wta_hay,nrow=1,ncol=3,common.legend = F,legend = "right")
ggsave('~/map_wta.jpeg',plot=map_wta_combined,width=16,height=8,dpi=300) # Fig 2

### Unit Cost ###
createMap_c <- function(.data, .name, .limits, .breaks, .value, maptitle){
  ggplot() +  
    geom_sf(data = lake_gsl, fill = "#c6dbef", color = "#9ecae1") +
    geom_tile(data = cdl_water_df, aes(x = x, y = y, fill = color), na.rm = TRUE, show.legend = FALSE) +
    scale_fill_manual(name = '', values = '#3182bd', labels = NULL, na.value = NA) +
    new_scale_fill() +  
    geom_tile(data = .data, aes(x = X, y = Y, fill = {{.value}})) +
    scale_fill_viridis_c(name = .name, limits = .limits, breaks = .breaks, guide = guide_colorbar(direction = "horizontal")) +
    geom_sf(data = basin_gsl, fill = NA, color = "gray30") +
    geom_sf_text(data = basin_gsl, aes(label = basin_name), color = "gray30") + 
    annotation_scale(location = "br", width_hint = 0.2, text_cex = 0.75, bar_cols = c("gray30","white")) +
    theme_minimal() +
    labs(title = maptitle, fill = "Value") +
    theme(
      plot.title = element_text(size = 16, hjust = 0.5, vjust = 1,face = "bold"),  # Adjust the size, face, and alignment
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_blank(),       # Remove x-axis text
      axis.text.y = element_blank(),       # Remove y-axis text
      axis.ticks = element_blank(),        # Remove axis ticks
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      panel.border = element_blank(),    
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 12),
      legend.direction = "horizontal",
      legend.key.height = unit(0.6, "cm"),    # Adjust height
      legend.key.width = unit(1.2, "cm")      # Adjust width
    ) 
}
createMap_c_noscale <- function(.data, .name, .limits, .breaks, .value, maptitle){
  ggplot() +  
    geom_sf(data = lake_gsl, fill = "#c6dbef", color = "#9ecae1") +
    geom_tile(data = cdl_water_df, aes(x = x, y = y, fill = color), na.rm = TRUE, show.legend = FALSE) +
    scale_fill_manual(name = '', values = '#3182bd', labels = NULL, na.value = NA) +
    new_scale_fill() +  
    geom_tile(data = .data, aes(x = X, y = Y, fill = {{.value}})) +
    scale_fill_viridis_c(name = .name, limits = .limits, breaks = .breaks,guide = guide_colorbar(direction = "horizontal")) +
    geom_sf(data = basin_gsl, fill = NA, color = "gray30") +
    geom_sf_text(data = basin_gsl, aes(label = basin_name), color = "gray30") + 
    theme_minimal() +
    labs(title = maptitle, fill = "Value") +
    theme(
      plot.title = element_text(size = 16, hjust = 0.5, vjust = 1,face = "bold"),  # Adjust the size, face, and alignment
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_blank(),       # Remove x-axis text
      axis.text.y = element_blank(),       # Remove y-axis text
      axis.ticks = element_blank(),        # Remove axis ticks
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      panel.border = element_blank(),    
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 12),
      legend.direction = "horizontal",
      legend.key.height = unit(0.6, "cm"),    # Adjust height
      legend.key.width = unit(1.2, "cm"),       # Adjust width
    ) 
}
cost_limits_fallow <- range(c(alfalfa_cost_fallow$value),na.rm = TRUE)
cost_breaks_fallow <- pretty(cost_limits_fallow)
map_cost_fallow <- createMap_c_noscale(.data = alfalfa_cost_fallow, maptitle = '(a) To Fallow',
                               .name = 'Cost ($/m3)', .limits = cost_limits_fallow, .breaks = cost_breaks_fallow, .value = value)
cost_limits_grain <- range(c(alfalfa_cost_grain$value),na.rm = TRUE)
cost_breaks_grain <- pretty(cost_limits_grain)
map_cost_grain <- createMap_c_noscale(.data = alfalfa_cost_grain, maptitle = '(b) To Spring Grains',
                                    .name = 'Cost ($/m3)', .limits = cost_limits_grain, .breaks = cost_breaks_grain, .value = value)
cost_limits_hay <- range(c(alfalfa_cost_hay$value),na.rm = TRUE)
cost_breaks_hay <- pretty(cost_limits_hay)
map_cost_hay <- createMap_c(.data = alfalfa_cost_hay, maptitle = '(c\u200B) To Other Hays', 
                          .name = 'Cost ($/m3)', .limits = cost_limits_hay, .breaks = cost_breaks_hay, .value = value)
map_cost_combined <- ggarrange(map_cost_fallow,map_cost_grain,map_cost_hay,nrow=1,ncol=3,common.legend = F,legend="bottom")
ggsave('~/map_cost.jpeg',plot=map_cost_combined,width=16,height=8,dpi=300)

map_final <- ggdraw() +
  draw_plot(map_cost_combined) +
  draw_plot(inset_map, x = 0.18, y = 0.13, width = 0.27, height = 0.27) +
  draw_plot(arrow_plot, x = 0.135, y = 0.16, width = 0.2, height = 0.2) 
ggsave('~/map_final1.jpeg',plot=map_final,width=16,height=8,dpi=300)

county_jordan <- st_intersection(subset(basin_gsl,basin_name=='JORDAN'), 
                                 subset(county,NAME=='JUAB'|NAME=='UTAH'|NAME=='WASATCH'))
cost_jordan <- alfalfa_wta %>% 
  filter(basin=='JORDAN RIVER', county=='JUAB'|county=='UTAH'|county=='WASATCH') %>% 
  rename(value = cost_Fallow) %>% 
  mutate(color = NA) %>% 
  mutate(bins=cut(value,quantile(cost_jordan$value,probs=seq(0,1,0.1)),include.lowest = TRUE)) %>% 
  dplyr::select(X,Y,value,color,bins)
#limits_jordan <- range(c(cost_jordan$value),na.rm = TRUE)
limits_jordan <- quantile(cost_jordan$value, probs = c(0, 0.99))
breaks_jordan <- seq(limits_jordan[1], limits_jordan[2], length.out = 5) %>% round(.,2)
inset_map <- ggplot() +
  geom_sf(data = county_jordan, fill = NA, color = "red", size = 1.5) +
  geom_sf_text(data = county_jordan, aes(label = NAME), color = "black", size = 3, fontface = "bold") + 
  geom_tile(data = cost_jordan, aes(x = X, y = Y, fill = value)) +
  scale_fill_gradientn(
    colors = c("#542788", "#4daf4a", '#ffff33','#e41a1c'),  # Darker, mid-range, and brighter colors
    name = '', 
    limits = limits_jordan, 
    breaks = breaks_jordan,
    guide = guide_colorbar(direction = "horizontal")
  ) +
  #scale_fill_viridis_c(name = '', limits = limits_jordan, breaks = breaks_jordan, option = 'inferno',guide = guide_colorbar(direction = "horizontal")) +
  theme_void() +
  theme(legend.position = c(0.55, -0.05),
        legend.justification = c("bottom"),
        legend.direction = "horizontal",
        legend.key.height = unit(0.2, "cm"),    # Adjust height
        legend.key.width = unit(0.5, "cm"),       # Adjust width
        legend.text = element_text(size = 7),   # Adjust text size
        )

arrow_plot <- ggplot() +
  annotate("segment", x = 0.13, xend = 0.1, y = 0.16, yend = 0.16, 
           arrow = arrow(length = unit(0.2, "inches"), type = "closed"), 
           color = "red", lwd = 1) + # Customize arrow properties here
  coord_fixed() +  # Use fixed coordinates to keep proportions
  theme_void()  

map_cost_fallow_inset <- ggdraw() +
  draw_plot(map_cost_fallow, hjust = 0, vjust = 0, scale = 1) +
  draw_plot(inset_map,halign = 0.96, valign = 1, scale = .33)
map_cost_combined <- ggarrange(map_cost_fallow_inset,map_cost_grain,map_cost_hay,nrow=1,ncol=3,common.legend = F, legend = "right")
ggsave('~/map_cost.jpeg',plot=map_cost_combined,width=16,height=8,dpi=300) # Fig 3

####### Plot Figs 5 and 6 #######
alfalfa_wta_pixel <- import('~/alfalfa_wta_pixel.csv')
alfalfa_wta_county <- import('~/alfalfa_wta_county.csv') 
alfalfa_wta_basin <- import('~/alfalfa_wta_basin.csv')
alfalfa_wta_county[is.na(alfalfa_wta_county)] <- 0
alfalfa_wta_basin[is.na(alfalfa_wta_basin)] <- 0

createMap_a <- function(.value, .data, maptitle){
  ggplot() +  
    geom_sf(data = lake_gsl, fill = "#c6dbef", color = "#9ecae1") +
    geom_tile(data = cdl_water_df, aes(x = x, y = y, fill = color), na.rm = TRUE, show.legend = FALSE) +
    scale_fill_manual(name = '', values = '#3182bd', labels = NULL, na.value = NA) +
    new_scale_fill() +  
    geom_tile(data = .data, aes(x = X, y = Y, fill = factor({{.value}}))) +
    scale_fill_manual(name = 'Projected Status', values = c("1" = "#1a9850","0" = "#d73027"), labels = c("1" = "Enrolled","0" = "Not enrolled")) +
    geom_sf(data = basin_gsl, fill = NA, color = "gray30") +
    geom_sf_text(data = basin_gsl, aes(label = basin_name), color = "gray30") + 
    annotation_scale(location = "br", width_hint = 0.2, text_cex = 0.75, bar_cols = c("gray30","white")) +
    theme_minimal() +
    labs(title = maptitle, fill = "Value") +
    theme(
      plot.title = element_markdown(size = 20, hjust = 0.5, vjust = 1,face = "bold"),  # Adjust the size, face, and alignment
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 14),
      legend.key.size = unit(1.2, 'cm')
    ) +
    coord_sf()
}
createMap_a_noscale <- function(.value, .data, maptitle){
  ggplot() +
    geom_tile(data = .data, aes(x = X, y = Y, fill = factor({{.value}}))) +
    scale_fill_manual(name = 'Projected Status', values = c("1" = "#1a9850","0" = "#d73027"), labels = c("1" = "Enrolled","0" = "Not enrolled")) +
    geom_sf(data = basin_gsl, fill = NA, color = "gray30") +
    geom_sf_text(data = basin_gsl, aes(label = basin_name), color = "gray30") + 
    #annotation_scale(location = "bl", width_hint = 0.2, text_cex = 0.75, bar_cols = c("gray30","white")) +
    theme_minimal() +
    labs(title = maptitle, fill = "Value") +
    theme(
      plot.title = element_markdown(size = 20, hjust = 0.5, vjust = 1,face = "bold"),  # Adjust the size, face, and alignment
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 14),
      legend.key.size = unit(1.2, 'cm')
    ) +
    coord_sf()
}
createMap_cutoff <- function(.data, .name, .limits, .breaks, .value, maptitle){
  ggplot() +  
    geom_sf(data = lake_gsl, fill = "#c6dbef", color = "#9ecae1") +
    geom_tile(data = cdl_water_df, aes(x = x, y = y, fill = color), na.rm = TRUE, show.legend = FALSE) +
    scale_fill_manual(name = '', values = '#3182bd', labels = NULL, na.value = NA) +
    new_scale_fill() +  
    geom_tile(data = .data, aes(x = X, y = Y, fill = {{.value}})) +
    scale_fill_viridis_c(name = .name, limits = .limits, breaks = .breaks) +
    geom_sf(data = basin_gsl, fill = NA, color = "gray30") +
    geom_sf_text(data = basin_gsl, aes(label = basin_name), color = "gray30") + 
    annotation_scale(location = "br", width_hint = 0.2, text_cex = 0.75, bar_cols = c("gray30","white")) +
    theme_minimal() +
    labs(title = maptitle, fill = "Value") +
    theme(
      plot.title = element_text(size = 22, hjust = 0.5, vjust = 1,face = "bold"),  # Adjust the size, face, and alignment
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 16),
      legend.key.size = unit(1.2, 'cm')
    ) +
    coord_sf()
}
createMap_cutoff_noscale <- function(.data, .name, .limits, .breaks, .value, maptitle){
  ggplot() +  
    geom_sf(data = lake_gsl, fill = "#c6dbef", color = "#9ecae1") +
    geom_tile(data = cdl_water_df, aes(x = x, y = y, fill = color), na.rm = TRUE, show.legend = FALSE) +
    scale_fill_manual(name = '', values = '#3182bd', labels = NULL, na.value = NA) +
    new_scale_fill() +  
    geom_tile(data = .data, aes(x = X, y = Y, fill = {{.value}})) +
    scale_fill_viridis_c(name = .name, limits = .limits, breaks = .breaks) +
    geom_sf(data = basin_gsl, fill = NA, color = "gray30") +
    geom_sf_text(data = basin_gsl, aes(label = basin_name), color = "gray30") + 
    theme_minimal() +
    labs(title = maptitle, fill = "Value") +
    theme(
      plot.title = element_text(size = 22, hjust = 0.5, vjust = 1,face = "bold"),  # Adjust the size, face, and alignment
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 16),
      legend.key.size = unit(1.2, 'cm')
    ) +
    coord_sf()
}

map_baseline <- createMap_a_noscale(maptitle = '(a) Projected Status', .value = enroll, .data = alfalfa_wta_pixel)
map_baseline_county <- createMap_a_noscale(maptitle = '(b) Lower 95% CI (County)', .value = enroll_upper, .data = alfalfa_wta_county)
map_baseline_basin <- createMap_a(maptitle = '(c\u200B) Lower 95% CI (Watershed)',.value = enroll_upper, .data = alfalfa_wta_basin)
map_combined_baseline <- ggarrange(map_baseline,map_baseline_county,map_baseline_basin,nrow = 1,ncol = 3,common.legend = T,legend="right")
ggsave('/Users/manli/Library/CloudStorage/Box-Box/Manuscripts/UT_WTA/Tex/map_cutoff_baseline.jpeg',plot=map_combined_baseline,
       width=16,height=8,dpi=300)

#map_baseline <- createMap_a_noscale(maptitle = '', .value = enroll, .data = alfalfa_wta_pixel)
#map_conserved <- createMap_a_noscale(maptitle = '', .value = enroll_a, .data = alfalfa_wta_pixel)
#map_combined <- ggarrange(map_baseline, map_conserved,nrow=1,ncol=2,common.legend=T,legend="right",labels=c('(a)','(b)'))

cutoff_limits <- range(c(alfalfa_wta_pixel$cutoff_wta),na.rm = TRUE)
cutoff_breaks <- pretty(cutoff_limits)
map_baseline_cutoff <- createMap_cutoff_noscale(.data = alfalfa_wta_pixel, maptitle = '(a) Baseline', .name = 'Payment ($/ha)', .limits = cutoff_limits, 
                                 .breaks = cutoff_breaks, .value = cutoff_wta)
cutoff_limits <- range(c(alfalfa_wta_pixel$cutoff_wta_upper),na.rm = TRUE)
cutoff_breaks <- pretty(cutoff_limits)
map_conserved_cutoff <- createMap_cutoff(.data = alfalfa_wta_pixel, maptitle = '(b) Conservative', .name = 'Payment ($/ha)', .limits = cutoff_limits, 
                                 .breaks = cutoff_breaks, .value = cutoff_wta_upper)
map_combined_cutoff <- ggarrange(map_baseline_cutoff, map_conserved_cutoff,nrow=1,ncol=2,common.legend = T,legend="right")
ggsave('~/map_cutoff_payment.jpeg',plot=map_combined_cutoff,width=16,height=10,dpi=300) # Fig 5  
 
# Baseline scenario: county and basin, and uniform
map_baseline <- createMap_a_noscale(maptitle = 'Projection (County, Watershed)', .value = enroll, .data = alfalfa_wta_county)
map_baseline_county <- createMap_a_noscale(maptitle = 'Lower Bound of 95% CI (County)', .value = enroll_upper, .data = alfalfa_wta_county)
map_baseline_basin <- createMap_a(maptitle = 'Lower Bound of 95% CI (Watershed)', .value = enroll_upper, .data = alfalfa_wta_basin)
map_combined_baseline <- ggarrange(map_baseline,map_baseline_county,map_baseline_basin,nrow = 1, ncol = 3,
          common.legend = T,legend="bottom",labels=c('(a)','(b)','(c)'))
ggsave('~/map_cutoff_baseline.jpeg',plot=map_combined_baseline,
       width=16,height=10,dpi=300) # Fig 6

# Conservative scenario: county, basin, and uniform
map_conserved_county <- createMap_a_noscale(maptitle = '(a) Projected Status (County)', .value = enroll_a, .data = alfalfa_wta_county)
map_conserved_basin <- createMap_a_noscale(maptitle = '(b) Projected Status (Watershed)', .value = enroll_a, .data = alfalfa_wta_basin)
map_conserved <- createMap_a(maptitle = '(c\u200B) Lower 95% CI (County, Watershed)', .value = enroll_a_upper, .data = alfalfa_wta_county)
map_combined_conserved <- ggarrange(map_conserved_county,map_conserved_basin,map_conserved,nrow = 1, ncol = 3,
                                    common.legend = T,legend="right")
ggsave('~/map_cutoff_conserved.jpeg',plot=map_combined_conserved,
       width=16,height=8,dpi=300) # Fig S3
