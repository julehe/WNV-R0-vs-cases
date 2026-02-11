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

df_wnnd <- read.csv("data/WNND_data_with_covariates.csv")

#-------------------------------------------------------------------------------
# Calculate percentage of cases occurring July-September 
# and temperature statistics for these months
#-------------------------------------------------------------------------------

df_wnnd$month <- month(df_wnnd$Date)
df_wnnd$year <- year(df_wnnd$Date)

# Calculate total observations
total_cases <- sum(df_wnnd$wnnd_sum)

# Count observations in July-September
july_sept_cases <- sum(
  (df_wnnd %>%
     filter(month %in% 7:9))$wnnd_sum
)

# Calculate percentage
percentage_july_sept_cases <- (july_sept_cases / total_cases) * 100

# Print result
percentage_july_sept_cases
total_cases - july_sept_cases

# temperature statistics
mean((df_wnnd %>%
  filter(month %in% 7:9))$temp_3mo_avg_pop_weight,
  na.rm = T)

df_wnnd_filtered <- df_wnnd %>%
                      group_by(NUTS_ID, year) %>%
                      filter(any(wnnd_sum>0))%>%
                      ungroup() %>%
                      filter(month %in% 7:9)

mean(df_wnnd_filtered$temp_3mo_avg_pop_weight, na.rm=T)
sum(df_wnnd_filtered$temp_3mo_avg_pop_weight>24.4, na.rm = T)/sum(!is.na(df_wnnd_filtered$temp_3mo_avg_pop_weight))

#-------------------------------------------------------------------------------
# Plot case counts by month and country
#-------------------------------------------------------------------------------

df_wnnd$country <- substr(df_wnnd$NUTS_ID, 1, 2)

country_counts <- df_wnnd %>% 
  group_by(country) %>%
  summarize(count = sum(wnnd_sum))

hyper_endemics <- country_counts %>% 
  filter(country %in% c("EL", "IT", "RO", "HU"))

sum(hyper_endemics$count)/sum(country_counts$count)

# Define a threshold for "Others"
threshold <- 100

# Replace low-count countries with "Others"
df_summary <- df_wnnd 
df_summary$country <- ifelse(df_summary$country %in% country_counts$country[country_counts$count >= threshold], 
                             df_summary$country, "XOthers")

new_labels <- c("EL" = "Greece", "ES" = "Spain", "HR" = "Croatia",
                "HU" = "Hungary", "IT" = "Italy", "RO" = "Romania",
                "XOthers" = "Others")

df_summary <- df_summary %>%
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
  theme_bw(base_size = 9) + 
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

df_cases_per_year <- df_wnnd %>% 
  group_by(year) %>%
  summarise(wnnd_sum = sum(wnnd_sum))

plot_cases_per_year <- ggplot(df_cases_per_year, aes(x = as.numeric(year), y = wnnd_sum)) +
  geom_col(fill = "steelblue") +
  geom_smooth(method = "lm", se = FALSE, color = "black", 
              linewidth = 0.8) +
  scale_x_continuous(breaks = as.numeric(df_cases_per_year$year)) +
  labs(
    title = "Total WNND cases per year",
    x = "Year",
    y = "Number of cases"
  ) +
  theme_bw(base_size = 9) +
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

# country borders of non-EU/EEA regions
EU_NUTS0 <- get_eurostat_geospatial(resolution = 10,
                                    nuts_level = 0, 
                                    year = 2021) 

EU_NUTS0 <- filter(EU_NUTS0, (CNTR_CODE %in% c("RS", "CH", "MK",
                                               "BA", "ME", "XK", "AL"))) 

# nuts3 regions
EU_NUTS3 <- st_read("data/NUTS/NUTS3_2024_with_UK_2021_res_10.shp")

# remove northern countries for better visual focus
EU_NUTS3 <- filter(EU_NUTS3, !(CNTR_CODE %in% c("UK", "NO", "SE", "IE",
                                                "EE", "LV", "LT", "FI", "IS"))) 

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
  theme_bw(base_size = 9) +
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
# Plot map with total number of cases per NUTS3 
#-------------------------------------------------------------------------------

df_wnnd_total  <- df_wnnd  %>%
  group_by(NUTS_ID) %>%
  summarise(wnnd_sum = sum(wnnd_sum))

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

EU_NUTS <- left_join(EU_NUTS, df_wnnd_total, by = "NUTS_ID")

EU_NUTS_zero <- EU_NUTS %>% filter(wnnd_sum == 0)
EU_NUTS_nonzero <- EU_NUTS %>% filter(wnnd_sum > 0) 

plot_total_cases <- ggplot() +
  geom_sf(data = EU_NUTS_zero, fill = "grey90", color = NA) +
  geom_sf(data = EU_NUTS_nonzero, aes(fill = wnnd_sum), color = NA) +
  geom_sf(data = EU_NUTS0_full, fill = NA, linewidth = 0.3) +
  scale_fill_gradientn(
    colours = c("#74add1", "#ffffbf", "#f46d43", "#a50026"),
    na.value = "white",
    name = "",
    guide = guide_colorbar(barwidth = 10, barheight = 0.5),
    transform = "log2"
  ) +
  labs(
    title = expression("Total cases of WNND per region (" * log[2]*"-scale)"),
    fill=""
  ) + 
  theme_bw(base_size = 9) +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box.margin = margin(t = 10),
    plot.margin = margin(t = 6, r = 6, b = 6, l = 20),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

plot_total_cases

#-------------------------------------------------------------------------------
# Save plot grid
#-------------------------------------------------------------------------------

plot_list = list(plot_cases_month_country, plot_nbr_years_map,
                 plot_cases_per_year, plot_total_cases)

plot_grid = cowplot::plot_grid(plotlist = plot_list, ncol=2, label_size = 12,
                               align = "h", axis = "b", labels = c('A', 'B',
                                                                   'C','D')) 

plot_grid

ggsave("Figures/intro_plots/Figure_1.tiff", 
       plot = plot_grid,
       width = 7.3, height = 7, 
       bg = "white",
       dpi = 300, units = "in", compression = "lzw")
ggsave("Figures/intro_plots/Figure_1.svg", 
       plot = plot_grid,
       width = 7.3, height = 7, 
       bg = "white",
       units = "in")

