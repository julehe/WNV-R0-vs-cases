# ==============================================================================
# Author: Julian Heidecke        
# Email: julian.heidecke@iwr.uni-heidelberg.de or julian.heidecke@gmail.com
# ==============================================================================

library(tidyverse)
library(ggplot2)
library(cowplot)
library(eurostat)
library(sf)

#-------------------------------------------------------------------------------
# Plot map with total incidence per NUTS3 
#-------------------------------------------------------------------------------

# Load monthly data and summarise over period

df_wnnd <- read.csv("data/WNND_data_with_covariates.csv")

df_wnnd_summary <- df_wnnd  %>%
  mutate(year = year(Date)) %>%
  group_by(NUTS_ID, year) %>%
  summarise(r0_from_seasonal_summer = mean(r0_from_seasonal_3mo_avg_pop_weight[7:9]),
            population = round(mean(population)),
            wnnd_sum = sum(wnnd_sum),
            .groups = "drop") %>%
  group_by(NUTS_ID) %>%
  summarise(r0_from_seasonal_summer = mean(r0_from_seasonal_summer),
            population = round(mean(population)),
            wnnd_sum = sum(wnnd_sum),
            .groups = "drop") %>%
  mutate(total_incidence_per_100K = 100000 * wnnd_sum / population)

# country borders of non-EU/EEA regions
EU_NUTS0 <- get_eurostat_geospatial(resolution = 10,
                                    nuts_level = 0, 
                                    year = 2021) 

EU_NUTS0 <- filter(EU_NUTS0, (CNTR_CODE %in% c("RS", "CH", "MK",
                                               "BA", "ME", "XK", "AL"))) 

# nuts3 regions
EU_NUTS3 <- st_read("data/NUTS/NUTS3_2024_with_UK_2021_res_10.shp")

# remove iceland for better visual focus
EU_NUTS3 <- filter(EU_NUTS3, !(CNTR_CODE %in% c("IS"))) 

EU_NUTS <- rbind(EU_NUTS0, EU_NUTS3)

plot(EU_NUTS$geometry)

# all country borders 
EU_NUTS0_full <- get_eurostat_geospatial(resolution = 10,
                                         nuts_level = 0, 
                                         year = 2021) 

EU_NUTS0_full <- filter(EU_NUTS0_full, CNTR_CODE %in% EU_NUTS$CNTR_CODE)

manual_bbox <- sf::st_bbox(c(
  xmin = -11,
  ymin = 33,
  xmax = 35,
  ymax = 70.5
), crs = st_crs(EU_NUTS))  

EU_NUTS0_full <- st_crop(EU_NUTS0_full, manual_bbox)

EU_NUTS <- left_join(EU_NUTS, df_wnnd_summary, by = "NUTS_ID")

EU_NUTS_zero <- EU_NUTS %>% filter(wnnd_sum == 0)
EU_NUTS_nonzero <- EU_NUTS %>% filter(wnnd_sum > 0) 

plot_total_incidence <- ggplot() +
  geom_sf(data = EU_NUTS_zero, fill = "grey90", color = NA) +
  geom_sf(data = EU_NUTS_nonzero, aes(fill = total_incidence_per_100K), color = NA) +
  geom_sf(data = EU_NUTS0_full, fill = NA, linewidth = 0.3) +
  scale_fill_gradientn(
    colours = c("#74add1", "#ffffbf", "#f46d43", "#a50026"),
    breaks = c(0.0625,0.5,4,32),
    labels = c("0.0625", "0.5", "4", "32"),
    na.value = "white",
    name = "",
    guide = guide_colorbar(
      barwidth = 0.5,
      barheight = 10,
      title.position = "top",
      title.hjust = 0.5
    ),
    transform = "log2"
  ) +
  labs(
    title = expression("Total incidence of WNND per 100K (" * log[2]*"-scale)"),
    fill=""
  ) + 
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.position = "right",
    legend.direction = "vertical",
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.background = element_rect(fill = "transparent", color = NA), 
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.box.background = element_rect(fill = "transparent", color = NA)
  )

plot_total_incidence

#-------------------------------------------------------------------------------
# Plot average suitability map 
#-------------------------------------------------------------------------------

# prepare data for counting number of years with wnnd presence in each region
df_nbr_years_wnnd <- df_wnnd %>%
  select(NUTS_ID, year, wnnd_sum) %>%
  group_by(NUTS_ID, year) %>%
  summarise(wnnd_pa = any(wnnd_sum>0),
            .groups = "drop") %>%
  group_by(NUTS_ID) %>%
  summarise(nbr_years_wnnd = sum(wnnd_pa)) %>%
  mutate(nbr_years_wnnd = as.character(ifelse(nbr_years_wnnd >6, 6, nbr_years_wnnd))) %>%
  mutate(nbr_years_wnnd = ifelse(nbr_years_wnnd == "6", "≥6",nbr_years_wnnd))

df_nbr_years_wnnd$nbr_years_wnnd <- factor(
  df_nbr_years_wnnd$nbr_years_wnnd,
  levels = c("0", "1", "2", "3", "4", "5", "≥6")
)

# define contiguous hotspot regions in greece, italy, romania, hungary based on the 
# number of years with WNND presence
nuts_hotspots <- (filter(df_nbr_years_wnnd, nbr_years_wnnd %in% c("≥6")) %>%
                    filter(!(NUTS_ID %in% c("ITG2G"))) %>%
                    filter(!(NUTS_ID %in% c("EL305"))) %>%
                    filter(!(NUTS_ID %in% c("RO211","RO213","RO125","RO126"))))$NUTS_ID

EU_NUTS_italy <- EU_NUTS3 %>%
  filter(NUTS_ID %in% grep("^IT", nuts_hotspots, value = TRUE)) %>%
  st_union() %>%
  st_sf(geometry = .) 

EU_NUTS_greece <- EU_NUTS3 %>%
  filter(NUTS_ID %in% grep("^EL", nuts_hotspots, value = TRUE)) %>%
  st_union() %>%
  st_sf(geometry = .)  

EU_NUTS_hungary <- EU_NUTS3 %>%
  filter(NUTS_ID %in% grep("^HU", nuts_hotspots, value = TRUE)) %>%
  st_union() %>%
  st_sf(geometry = .)  

EU_NUTS_romania <- EU_NUTS3 %>%
  filter(NUTS_ID %in% grep("^RO", nuts_hotspots, value = TRUE)) %>%
  st_union() %>%
  st_sf(geometry = .)  

EU_NUTS_highlight <- rbind(EU_NUTS_italy, EU_NUTS_greece, 
                           EU_NUTS_hungary, EU_NUTS_romania)

R0_map <- ggplot() +
  geom_sf(data = EU_NUTS, aes(fill = r0_from_seasonal_summer), color = NA) +
  geom_sf(data = EU_NUTS0_full, fill = NA, linewidth = 0.3) +
  geom_sf(data = EU_NUTS_highlight, fill = NA, color = "black", linewidth = 0.4) +
  scale_fill_gradientn(
    colours = c("#74add1", "#ffffbf", "#f46d43", "#a50026"),
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    na.value = "white",
    name = "",
    guide = guide_colorbar(
      barwidth = 0.5,
      barheight = 10,
      title.position = "top",
      title.hjust = 0.5
    )
  ) +
  labs(
    title = expression(
      "Average Jul-Sep " * R[0]^{"rel"} * " in 2010–2023"
    ),
    fill=""
  ) + 
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.position = "right",
    legend.direction = "vertical",
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.background = element_rect(fill = "transparent", color = NA), 
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.box.background = element_rect(fill = "transparent", color = NA)
  )

R0_map

#-------------------------------------------------------------------------------
# Save maps
#-------------------------------------------------------------------------------

plot_total_incidence

ggsave("Figures/spatial_temporal_comparison/Incidence_map.tiff", 
       plot = plot_total_incidence,
       width = 4, height = 5, 
       dpi = 300, units = "in", compression = "lzw")
ggsave("Figures/spatial_temporal_comparison/Incidence_map.svg", 
       plot = plot_total_incidence,
       width = 4, height = 5, 
       units = "in")

R0_map

ggsave("Figures/spatial_temporal_comparison/R0_map.tiff", 
       plot = R0_map,
       width = 4, height = 5, 
       dpi = 300, units = "in", compression = "lzw")
ggsave("Figures/spatial_temporal_comparison/R0_map.svg", 
       plot = R0_map,
       width = 4, height = 5, 
       units = "in")


#-------------------------------------------------------------------------------
# Time series of monthly cases in hotspots
#-------------------------------------------------------------------------------

df_wnnd <- df_wnnd %>% dplyr::select(NUTS_ID, Date, wnnd_sum, incidence_per_100K,
                                     temp_3mo_avg_pop_weight, 
                                     r0_from_seasonal_3mo_avg_pop_weight, 
                                     population)

### Romania

data_RO <- df_wnnd %>% filter(NUTS_ID %in% grep("^RO", nuts_hotspots, value = TRUE))

data_RO$Date <- as.Date(data_RO$Date)

data_RO <- data_RO %>%
  group_by(Date) %>%
  summarise(wnnd_sum = sum(wnnd_sum,na.rm=T),
            r0_from_seasonal_3mo_avg = sum(r0_from_seasonal_3mo_avg_pop_weight * population, na.rm = TRUE) / 
              sum(population[!is.na(r0_from_seasonal_3mo_avg_pop_weight)], na.rm = TRUE)
  )

data_RO$Date <- as.Date(data_RO$Date)

data_RO$Year <- format(data_RO$Date, "%Y")

data_RO$wnnd_sum <- data_RO$wnnd_sum/max(data_RO$wnnd_sum)

plot_romania<-ggplot(data_RO, aes(x = Date, y = wnnd_sum)) +
  geom_col(fill = "gray85", color="black",linewidth = 0.15,width=30) +
  geom_line(aes(y = r0_from_seasonal_3mo_avg), color = "#a50026", linewidth = 0.35) +
  geom_step(color = "black", direction = "mid", linewidth = 0.25) +
  scale_x_date(
    breaks = as.Date(c("2010-01-01", "2014-01-01", "2018-01-01", "2022-01-01")),
    labels = c("2010", "2014", "2018", "2022"),
    expand = c(0.01, 0.01)
  ) +
  labs(
    title = "Southern Romania",
    x = "",
    y = expression("Std. case count and " * R[0]^{"rel"})
  ) +
  theme_classic(base_size = 7) + 
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_rect(fill = "transparent", color = NA), 
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.box.background = element_rect(fill = "transparent", color = NA)
  )

plot_romania

### Hungary

data_HU <- df_wnnd %>% filter(NUTS_ID %in% grep("^HU", nuts_hotspots, value = TRUE))

data_HU$Date <- as.Date(data_HU$Date)

data_HU <- data_HU %>%
  group_by(Date) %>%
  summarise(wnnd_sum = sum(wnnd_sum,na.rm=T),
            r0_from_seasonal_3mo_avg = sum(r0_from_seasonal_3mo_avg_pop_weight * population, na.rm = TRUE) / 
              sum(population[!is.na(r0_from_seasonal_3mo_avg_pop_weight)], na.rm = TRUE)
  )

data_HU$Date <- as.Date(data_HU$Date)

data_HU$Year <- format(data_HU$Date, "%Y")

data_HU$wnnd_sum <- data_HU$wnnd_sum/max(data_HU$wnnd_sum)

plot_hungary<-ggplot(data_HU, aes(x = Date, y = wnnd_sum)) +
  geom_col(fill = "gray85", color="black",linewidth = 0.15,width=30) +
  geom_line(aes(y = r0_from_seasonal_3mo_avg), color = "#a50026", linewidth = 0.35) +
  geom_step(color = "black", direction = "mid", linewidth = 0.25) +
  scale_x_date(
    breaks = as.Date(c("2010-01-01", "2014-01-01", "2018-01-01", "2022-01-01")),
    labels = c("2010", "2014", "2018", "2022"),
    expand = c(0.01, 0.01)
  ) +
  labs(
    title = "Eastern Hungary",
    x = "",
    y = expression("Std. case count and " * R[0]^{"rel"})
  ) +
  theme_classic(base_size = 7) + 
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_rect(fill = "transparent", color = NA), 
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.box.background = element_rect(fill = "transparent", color = NA)
  )

plot_hungary

### Italy

data_NIT <- df_wnnd %>% filter(NUTS_ID %in% grep("^IT", nuts_hotspots, value = TRUE))

data_NIT$Date <- as.Date(data_NIT$Date)

data_NIT <- data_NIT %>%
  group_by(Date) %>%
  summarise(wnnd_sum = sum(wnnd_sum,na.rm=T),
            r0_from_seasonal_3mo_avg = sum(r0_from_seasonal_3mo_avg_pop_weight * population, na.rm = TRUE) / 
              sum(population[!is.na(r0_from_seasonal_3mo_avg_pop_weight)], na.rm = TRUE)
  )

data_NIT$Date <- as.Date(data_NIT$Date)

data_NIT$Year <- format(data_NIT$Date, "%Y")

data_NIT$wnnd_sum <- data_NIT$wnnd_sum/max(data_NIT$wnnd_sum)

plot_italy<-ggplot(data_NIT, aes(x = Date, y = wnnd_sum)) +
  geom_col(fill = "gray85", color="black",linewidth = 0.15,width=30) +
  geom_line(aes(y = r0_from_seasonal_3mo_avg), color = "#a50026", linewidth = 0.35) + 
  geom_step(color = "black", direction = "mid", linewidth = 0.25) +
  scale_x_date(
    breaks = as.Date(c("2010-01-01", "2014-01-01", "2018-01-01", "2022-01-01")),
    labels = c("2010", "2014", "2018", "2022"),
    expand = c(0.01, 0.01)
  ) +
  labs(
    title = "Northern Italy",
    x = "",
    y = expression("Std. case count and " * R[0]^{"rel"})
  ) +
  theme_classic(base_size = 7) + 
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_rect(fill = "transparent", color = NA), 
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.box.background = element_rect(fill = "transparent", color = NA)
  )

plot_italy

### Greece

data_main_greece <- df_wnnd %>% filter(NUTS_ID %in% grep("^EL", nuts_hotspots, value = TRUE))

data_main_greece$Date <- as.Date(data_main_greece$Date)

data_main_greece <- data_main_greece %>%
  group_by(Date) %>%
  summarise(wnnd_sum = sum(wnnd_sum,na.rm=T),
            r0_from_seasonal_3mo_avg = sum(r0_from_seasonal_3mo_avg_pop_weight * population, na.rm = TRUE) / 
              sum(population[!is.na(r0_from_seasonal_3mo_avg_pop_weight)], na.rm = TRUE)
  )

data_main_greece$Date <- as.Date(data_main_greece$Date)
data_main_greece$Year <- format(data_main_greece$Date, "%Y")

data_main_greece$wnnd_sum <- data_main_greece$wnnd_sum/max(data_main_greece$wnnd_sum)

plot_greece<-ggplot(data_main_greece, aes(x = Date, y = wnnd_sum)) +
  geom_col(fill = "gray85", color="black",linewidth = 0.15,width=30) +
  geom_line(aes(y = r0_from_seasonal_3mo_avg), color = "#a50026", linewidth = 0.35) + 
  geom_step(color = "black", direction = "mid", linewidth = 0.25) +
  scale_x_date(
    breaks = as.Date(c("2010-01-01", "2014-01-01", "2018-01-01", "2022-01-01")),
    labels = c("2010", "2014", "2018", "2022"),
    expand = c(0.01, 0.01)
  ) +
  labs(
    title = "Northern Greece",
    x = "",
    y = expression("Std. case count and " * R[0]^{"rel"})
  ) +
  theme_classic(base_size = 7) + 
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_rect(fill = "transparent", color = NA), 
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.box.background = element_rect(fill = "transparent", color = NA)
  )

plot_greece

#-------------------------------------------------------------------------------
# Save plots
#-------------------------------------------------------------------------------

ggsave("Figures/spatial_temporal_comparison/romania.tiff", 
       plot = plot_romania,
       width = 2.8, height = 1.5, 
       dpi = 300, units = "in", compression = "lzw")
ggsave("Figures/spatial_temporal_comparison/romania.svg", 
       plot = plot_romania,
       width = 2.8, height = 1.5, 
       units = "in")

ggsave("Figures/spatial_temporal_comparison/hungary.tiff", 
       plot = plot_hungary,
       width = 2.8, height = 1.5, 
       dpi = 300, units = "in", compression = "lzw")
ggsave("Figures/spatial_temporal_comparison/hungary.svg", 
       plot = plot_hungary,
       width = 2.8, height = 1.5, 
       units = "in")

ggsave("Figures/spatial_temporal_comparison/italy.tiff", 
       plot = plot_italy,
       width = 2.8, height = 1.5, 
       dpi = 300, units = "in", compression = "lzw")
ggsave("Figures/spatial_temporal_comparison/italy.svg", 
       plot = plot_italy,
       width = 2.8, height = 1.5, 
       units = "in")

ggsave("Figures/spatial_temporal_comparison/greece.tiff", 
       plot = plot_greece,
       width = 2.8, height = 1.5, 
       dpi = 300, units = "in", compression = "lzw")
ggsave("Figures/spatial_temporal_comparison/greece.svg", 
       plot = plot_greece,
       width = 2.8, height = 1.5, 
       units = "in")
