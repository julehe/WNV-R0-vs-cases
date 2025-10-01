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
# Load monthly data 
#-------------------------------------------------------------------------------

# The covariate data is openly available from:
# - Copernicus (https://cds.climate.copernicus.eu/datasets/derived-era5-land-daily-statistics?tab=overview)
# - Eurostat (https://ec.europa.eu/eurostat/de/)
# - CORINE (https://land.copernicus.eu/en/products/corine-land-cover)
# Restrictions apply to the Human WNND data which is only available upon request
# at TESSy (https://www.ecdc.europa.eu/en/publications-data/european-surveillance-system-tessy)
df_wnnd_ts <- read.csv("....csv")

#-------------------------------------------------------------------------------
# Calculate percentage of cases occurring July-September
#-------------------------------------------------------------------------------

df_wnnd_ts$month <- month(df_wnnd_ts$Date)
df_wnnd_ts$year <- year(df_wnnd_ts$Date)

# Calculate total observations
total_cases <- sum(df_wnnd_ts$wnnd_sum)

# Count observations in July-September
july_sept_cases <- sum(
  (df_wnnd_ts %>%
     filter(month %in% 7:9))$wnnd_sum
)

# Calculate percentage
percentage_july_sept_cases <- (july_sept_cases / total_cases) * 100

# Print result
percentage_july_sept_cases
total_cases - july_sept_cases

#-------------------------------------------------------------------------------
# Plot case counts by month and country
#-------------------------------------------------------------------------------

df_wnnd_ts$month <- month(df_wnnd_ts$Date)
df_wnnd_ts$country <- substr(df_wnnd_ts$NUTS_ID, 1, 2)

country_counts <- df_wnnd_ts %>% 
  group_by(country) %>%
  summarize(count = sum(wnnd_sum))

hyper_endemics <- country_counts %>% 
  filter(country %in% c("EL", "IT", "RO", "HU"))

sum(hyper_endemics$count)/sum(country_counts$count)

# Define a threshold for "Others"
threshold <- 100

# Replace low-count countries with "Others"
df_wnnd_ts$country <- ifelse(df_wnnd_ts$country %in% country_counts$country[country_counts$count >= threshold], 
                             df_wnnd_ts$country, "XOthers")

new_labels <- c("EL" = "Greece", "ES" = "Spain", "HR" = "Croatia",
                "HU" = "Hungary", "IT" = "Italy", "RO" = "Romania",
                "XOthers" = "Others")

df_summary <- df_wnnd_ts %>%
  group_by(month, country) %>%
  summarise(total_cases = sum(wnnd_sum), .groups = "drop")

df_summary <- df_summary %>%
  uncount(total_cases)  

plot_cases_month_country <- ggplot(df_summary, aes(x = month, fill = country)) +
  geom_histogram(binwidth = 1, color = "black", boundary = 0.5, position = "stack") +
  scale_x_continuous(breaks = 1:12, labels = month.abb, name = "Month") +
  labs(
    title = "WNND cases by month and country", 
    fill = "Country") +
  scale_fill_brewer(palette = "Set3", 
                    labels = new_labels,
                    guide = guide_legend(nrow=2)) +  
  theme_bw(base_size = 10) + 
  theme(panel.grid = element_blank(),
        legend.position = "bottom",  
        legend.direction = "horizontal",
        axis.title.x = element_blank(),
        axis.title.y = element_blank()
  )

plot_cases_month_country

#-------------------------------------------------------------------------------
# Plot total number of cases per year
#-------------------------------------------------------------------------------

df_cases_per_year <- df_wnnd_ts %>% 
  group_by(year) %>%
  summarise(wnnd_sum = sum(wnnd_sum))

plot_cases_per_year <- ggplot(df_cases_per_year, aes(x = as.numeric(year), y = wnnd_sum)) +
  geom_col(fill = "steelblue") +
  geom_smooth(method = "lm", se = FALSE, color = "black", #lty="dashed",
              linewidth = 0.8) +
  scale_x_continuous(breaks = as.numeric(df_cases_per_year$year)) +
  labs(
    title = "Total WNND cases per year",
    x = "Year",
    y = "Number of cases"
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )

plot_cases_per_year

#-------------------------------------------------------------------------------
# Plot map showing nbr of years with WNND presence for each NUTS3 region
#-------------------------------------------------------------------------------

EU_NUTS0 <- get_eurostat_geospatial(resolution = 10,
                                    nuts_level = 0, 
                                    year = 2021) 

EU_NUTS0 <- filter(EU_NUTS0, (CNTR_CODE %in% c("RS", "CH", "MK",
                                               "BA", "ME", "XK", "AL"))) 

EU_NUTS3 <- st_read("data/NUTS/NUTS3_2024_with_UK_2021_res_10.shp")

EU_NUTS3 <- filter(EU_NUTS3, !(CNTR_CODE %in% c("UK", "NO", "SE", "IE",
                                                "EE", "LV", "LT", "FI", "IS"))) 

EU_NUTS <- rbind(EU_NUTS0, EU_NUTS3)

plot(EU_NUTS$geometry)

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


#-------------------------------------------------------------------------------
# Load yearly data 
#-------------------------------------------------------------------------------

# The covariate data is openly available from:
# - Copernicus (https://cds.climate.copernicus.eu/datasets/derived-era5-land-daily-statistics?tab=overview)
# - Eurostat (https://ec.europa.eu/eurostat/de/)
# - CORINE (https://land.copernicus.eu/en/products/corine-land-cover)
# Restrictions apply to the Human WNND data which is only available upon request
# at TESSy (https://www.ecdc.europa.eu/en/publications-data/european-surveillance-system-tessy)
df_nbr_years_wnnd <- read.csv("....csv")

df_nbr_years_wnnd <- df_nbr_years_wnnd %>%
  filter(Year > 2009) %>%
  group_by(NUTS_ID) %>%
  summarise(nbr_years_wnnd = sum(wnnd_pa)) %>%
  mutate(nbr_years_wnnd = as.character(ifelse(nbr_years_wnnd >6, 6, nbr_years_wnnd))) %>%
  mutate(nbr_years_wnnd = ifelse(nbr_years_wnnd == "6", "≥6",nbr_years_wnnd))

df_nbr_years_wnnd$nbr_years_wnnd <- factor(
  df_nbr_years_wnnd$nbr_years_wnnd,
  levels = c("0", "1", "2", "3", "4", "5", "≥6")
)

EU_NUTS <- left_join(EU_NUTS, df_nbr_years_wnnd, by = "NUTS_ID")

EU_NUTS$nbr_years_wnnd <- factor(EU_NUTS$nbr_years_wnnd, levels = c("0", "1", "2", "3", "4", "5", "≥6"))

# Define colors: grey for "0", viridis for the rest
levels_all <- sort(unique(na.omit(EU_NUTS$nbr_years_wnnd)))
wnnd_levels <- setdiff(levels_all, "0")

colors <- c(
  "0" = "grey90",  
  setNames(viridis::viridis(length(wnnd_levels), option = "D"), wnnd_levels)
)

# Plot
plot_nbr_years_map <- ggplot(EU_NUTS) +
  geom_sf(aes(fill = nbr_years_wnnd), color = NA) +
  geom_sf(data = EU_NUTS0_full, fill = NA, linewidth = 0.3) +
  scale_fill_manual(
    values = colors,
    na.value = "white",
    drop = FALSE,
    guide = guide_legend(override.aes = list(colour = "black"))
  ) +
  labs(
    title = "Number of years with WNND cases",
    fill=""
  ) + 
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom",  
    legend.direction = "horizontal",
    legend.box.margin = margin(t = 10, r = 0, b = 0, l = 0), # 0 is default
    plot.margin = margin(t = 6, r = 6, b = 6, l = 20), # 6 is default
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank()
  )

plot_nbr_years_map

#-------------------------------------------------------------------------------
# Plot map with total incidence and total number of cases per NUTS3 
#-------------------------------------------------------------------------------

#-------------------
# Load yearly data 
#-------------------

# The covariate data is openly available from:
# - Copernicus (https://cds.climate.copernicus.eu/datasets/derived-era5-land-daily-statistics?tab=overview)
# - Eurostat (https://ec.europa.eu/eurostat/de/)
# - CORINE (https://land.copernicus.eu/en/products/corine-land-cover)
# Restrictions apply to the Human WNND data which is only available upon request
# at TESSy (https://www.ecdc.europa.eu/en/publications-data/european-surveillance-system-tessy)
df_wnnd <- read.csv("....csv")

df_wnnd_avg  <- df_wnnd  %>%
  group_by(NUTS_ID) %>%
  summarise(temp_summer=mean(temp_summer),
            wnnd_sum = sum(wnnd_sum),
            r0_post_from_seasonal_summer = mean(r0_post_from_seasonal_summer),
            incidence_per_100K = mean(incidence_per_100K),
            population= round(mean(population)))

EU_NUTS0 <- get_eurostat_geospatial(resolution = 10,
                                    nuts_level = 0, 
                                    year = 2021) 

EU_NUTS0 <- filter(EU_NUTS0, (CNTR_CODE %in% c("RS", "CH", "MK",
                                               "BA", "ME", "XK", "AL"))) 

EU_NUTS3 <- st_read("data/NUTS/NUTS3_2024_with_UK_2021_res_10.shp")

EU_NUTS3 <- filter(EU_NUTS3, !(CNTR_CODE %in% c("UK", "NO", "SE", "IE",
                                                "EE", "LV", "LT", "FI", "IS"))) 

EU_NUTS <- rbind(EU_NUTS0, EU_NUTS3)

plot(EU_NUTS$geometry)

EU_NUTS0_full <- get_eurostat_geospatial(resolution = 10,
                                         nuts_level = 0, 
                                         year = 2021) 

EU_NUTS0_full <- filter(EU_NUTS0_full, CNTR_CODE %in% EU_NUTS$CNTR_CODE)
EU_NUTS0_full <- st_crop(EU_NUTS0_full, manual_bbox)

EU_NUTS <- left_join(EU_NUTS, df_wnnd_avg, by = "NUTS_ID")

EU_NUTS_zero <- EU_NUTS %>% filter(wnnd_sum == 0)
EU_NUTS_nonzero <- EU_NUTS %>% filter(wnnd_sum > 0) %>%
  mutate(total_incidence = 100000 * wnnd_sum / population)

plot_total_incidence<- ggplot() +
  geom_sf(data = EU_NUTS_zero, fill = "grey90", color = NA) +
  geom_sf(data = EU_NUTS_nonzero, aes(fill = total_incidence), color = NA) +
  geom_sf(data = EU_NUTS0_full, fill = NA, linewidth = 0.3) +
  scale_fill_gradientn(
    colours = c("#74add1", "#ffffbf", "#f46d43", "#a50026"),
    na.value = "white",
    name = "",
    guide = guide_colorbar(barwidth = 10, barheight = 0.5)
  ) +
  labs(
    title = "Total incidence of WNND per 100K",
    fill=""
  ) + 
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box.margin = margin(t = 10),
    plot.margin = margin(t = 6, r = 6, b = 6, l = 20),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

plot_total_cases<-ggplot() +
  geom_sf(data = EU_NUTS_zero, fill = "grey90", color = NA) +
  geom_sf(data = EU_NUTS_nonzero, aes(fill = wnnd_sum), color = NA) +
  geom_sf(data = EU_NUTS0_full, fill = NA, linewidth = 0.3) +
  scale_fill_gradientn(
    colours = c("#74add1", "#ffffbf", "#f46d43", "#a50026"),
    breaks=c(50,100,150,200),
    na.value = "white",
    name = "",
    guide = guide_colorbar(barwidth = 10, barheight = 0.5)
  ) +
  labs(
    title = "Total cases of WNND per region",
    fill=""
  ) + 
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box.margin = margin(t = 10),
    plot.margin = margin(t = 6, r = 6, b = 6, l = 20),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

#-------------------------------------------------------------------------------
# Save plot grid
#-------------------------------------------------------------------------------

plot_list = list(plot_cases_month_country, plot_nbr_years_map,
                 plot_cases_per_year ,plot_total_cases)

plot_grid = cowplot::plot_grid(plotlist = plot_list, ncol=2, label_size = 12,
                               align = "h", axis = "b", labels = c('A', 'B',
                                                                   'C','D')) 

plot_grid

ggsave(".../main_intro.tiff", 
       plot = plot_grid,
       width = 8.3, height = 8, 
       bg = "white",
       dpi = 600, units = "in", compression = "lzw")

plot_total_incidence

ggsave(".../Incidence_map.tiff", 
       plot = plot_total_incidence,
       width = 4.15, height = 4, 
       dpi = 600, units = "in", compression = "lzw")

#-------------------------------------------------------------------------------
# Plot average suitability map 
#-------------------------------------------------------------------------------

EU_NUTS0 <- get_eurostat_geospatial(resolution = 10,
                                    nuts_level = 0, 
                                    year = 2021) 

EU_NUTS0 <- filter(EU_NUTS0, (CNTR_CODE %in% c("RS", "CH", "MK",
                                               "BA", "ME", "XK", "AL"))) 

EU_NUTS3 <- st_read("data/NUTS/NUTS3_2024_with_UK_2021_res_10.shp")

EU_NUTS <- rbind(EU_NUTS0, EU_NUTS3)

EU_NUTS <- left_join(EU_NUTS, df_wnnd_avg, by = "NUTS_ID")

EU_NUTS1 <- get_eurostat_geospatial(resolution = 10,
                                    nuts_level = 1, 
                                    year = 2024) 

# some hotspot nuts3 in greece, italy, romania, hungary
nuts_hotspots <- (filter(df_nbr_years_wnnd, nbr_years_wnnd %in% c("≥6")) %>%
                    filter(!(NUTS_ID %in% c("ITG2G"))) %>%
                    filter(!(NUTS_ID %in% c("EL612","EL305"))) %>%
                    filter(!(NUTS_ID %in% c("RO211","RO213","RO125","RO126"))))$NUTS_ID

EU_NUTS_nit <- EU_NUTS3 %>%
  filter(NUTS_ID %in% grep("^IT", nuts_hotspots, value = TRUE)) %>%
  st_union() %>%
  st_sf(geometry = .) 

EU_NUTS_maingreece <- EU_NUTS3 %>%
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

EU_NUTS_highlight <- rbind(EU_NUTS_nit, EU_NUTS_maingreece, 
                           EU_NUTS_hungary, EU_NUTS_romania)

EU_NUTS0_full <- get_eurostat_geospatial(resolution = 10,
                                         nuts_level = 0, 
                                         year = 2021) 

EU_NUTS0_full <- filter(EU_NUTS0_full, CNTR_CODE %in% EU_NUTS$CNTR_CODE)
EU_NUTS0_full <- st_crop(EU_NUTS0_full, manual_bbox)

R0_map<-ggplot() +
  geom_sf(data = EU_NUTS, aes(fill = r0_post_from_seasonal_summer), color = NA) +
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
    title = expression("Average summer " * R[0]^{"rel"} * " in 2010–2023"),
    fill=""
  ) + 
  theme_bw(base_size = 8) +
  theme(
    plot.title = element_text(hjust = 0.5), 
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.position = "right",
    legend.direction = "vertical",
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

R0_map

ggsave(".../R0_map_new.tiff", 
       plot = R0_map,
       width = 4, height = 5, 
       dpi = 600, units = "in", compression = "lzw")

#-------------------------------------------------------------------------------
# Time series of monthly cases in hotspots
#-------------------------------------------------------------------------------

#---------------------
# Load monthly data 
#---------------------

# The covariate data is openly available from:
# - Copernicus (https://cds.climate.copernicus.eu/datasets/derived-era5-land-daily-statistics?tab=overview)
# - Eurostat (https://ec.europa.eu/eurostat/de/)
# - CORINE (https://land.copernicus.eu/en/products/corine-land-cover)
# Restrictions apply to the Human WNND data which is only available upon request
# at TESSy (https://www.ecdc.europa.eu/en/publications-data/european-surveillance-system-tessy)
df_wnnd <- read.csv(".....csv")

df_wnnd <- df_wnnd %>% dplyr::select(NUTS_ID, Date, wnnd_sum, incidence_per_100K,
                                     temp_3mo_avg, r0_from_seasonal_3mo_avg, area_km2)

### Romania

data_RO <- df_wnnd %>% filter(NUTS_ID %in% grep("^RO", nuts_hotspots, value = TRUE))

data_RO$Date <- as.Date(data_RO$Date)

data_RO <- data_RO %>%
  group_by(Date) %>%
  summarise(wnnd_sum = sum(wnnd_sum,na.rm=T),
            r0_from_seasonal_3mo_avg = sum(r0_from_seasonal_3mo_avg * area_km2, na.rm = TRUE) / 
              sum(area_km2[!is.na(r0_from_seasonal_3mo_avg)], na.rm = TRUE)
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

ggsave(".../romania_new.tiff", 
       plot = plot_romania,
       width = 2.8, height = 1.5, 
       dpi = 600, units = "in", compression = "lzw")

### Hungary

data_HU <- df_wnnd %>% filter(NUTS_ID %in% grep("^HU", nuts_hotspots, value = TRUE))

data_HU$Date <- as.Date(data_HU$Date)

data_HU <- data_HU %>%
  group_by(Date) %>%
  summarise(wnnd_sum = sum(wnnd_sum,na.rm=T),
            r0_from_seasonal_3mo_avg = sum(r0_from_seasonal_3mo_avg * area_km2, na.rm = TRUE) / 
              sum(area_km2[!is.na(r0_from_seasonal_3mo_avg)], na.rm = TRUE)
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

ggsave(".../hungary_new.tiff", 
       plot = plot_hungary,
       width = 2.8, height = 1.5, 
       dpi = 600, units = "in", compression = "lzw")

### Italy

data_NIT <- df_wnnd %>% filter(NUTS_ID %in% grep("^IT", nuts_hotspots, value = TRUE))

data_NIT$Date <- as.Date(data_NIT$Date)

data_NIT <- data_NIT %>%
  group_by(Date) %>%
  summarise(wnnd_sum = sum(wnnd_sum,na.rm=T),
            r0_from_seasonal_3mo_avg = sum(r0_from_seasonal_3mo_avg * area_km2, na.rm = TRUE) / 
              sum(area_km2[!is.na(r0_from_seasonal_3mo_avg)], na.rm = TRUE)
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

ggsave(".../italy_new.tiff", 
       plot = plot_italy,
       width = 2.8, height = 1.5, 
       dpi = 600, units = "in", compression = "lzw")

### Greece

data_main_greece <- df_wnnd %>% filter(NUTS_ID %in% grep("^EL", nuts_hotspots, value = TRUE))

data_main_greece$Date <- as.Date(data_main_greece$Date)

data_main_greece <- data_main_greece %>%
  group_by(Date) %>%
  summarise(wnnd_sum = sum(wnnd_sum,na.rm=T),
            r0_from_seasonal_3mo_avg = sum(r0_from_seasonal_3mo_avg * area_km2, na.rm = TRUE) / 
              sum(area_km2[!is.na(r0_from_seasonal_3mo_avg)], na.rm = TRUE)
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

ggsave(".../greece_new.tiff", 
       plot = plot_greece,
       width = 2.8, height = 1.5, 
       dpi = 600, units = "in", compression = "lzw")

#-------------------------------------------------------------------------------
# Plot suitability map for each year
#-------------------------------------------------------------------------------

# The covariate data is openly available from:
# - Copernicus (https://cds.climate.copernicus.eu/datasets/derived-era5-land-daily-statistics?tab=overview)
# - Eurostat (https://ec.europa.eu/eurostat/de/)
# - CORINE (https://land.copernicus.eu/en/products/corine-land-cover)
# Restrictions apply to the Human WNND data which is only available upon request
# at TESSy (https://www.ecdc.europa.eu/en/publications-data/european-surveillance-system-tessy)
df_wnnd <- read.csv("....csv")

for(year in 2010:2023){
  df_temp <- filter(df_wnnd, Year == year)
  
  NUTS_temp <- rbind(EU_NUTS0, EU_NUTS3)
  
  NUTS_temp <- left_join(NUTS_temp, df_temp, by = "NUTS_ID")
  
  NUTS_cases <- filter(NUTS_temp, wnnd_pa == 1)
  NUTS_cases_centroids <- st_centroid(NUTS_cases)
  
  map_temp <- ggplot() +
    geom_sf(data = NUTS_temp, aes(fill = r0_post_from_seasonal_summer), color = NA) +
    geom_sf(data = EU_NUTS0_full, fill = NA, linewidth = 0.3) +
    geom_sf(data = NUTS_cases_centroids, color = "black", size = 1.5, shape = 21, fill = "yellow") +
    scale_fill_gradientn(
      colours = c("#74add1", "#ffffbf", "#f46d43", "#a50026"),
      limits = c(0, 1),
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      na.value = "white",
      name = bquote("Summer " * R[0]^{"rel"} * " in " * .(year)),
      guide = guide_colorbar(
        barwidth = 10,
        barheight = 0.5,
        title.position = "top",
        title.hjust = 0.5
      )
    ) +
    labs(
      title = "",
      fill=""
    ) + 
    theme_bw(base_size = 9) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      legend.position = "bottom",
      legend.direction = "horizontal",
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )
  
  print(map_temp)
  
}
