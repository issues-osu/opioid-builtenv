options(scipen = 999)  
options(tigris_use_cache = TRUE)

library(tidycensus)
library(tigris)
library(dplyr)
library(sf)
library(terra)
library(readr)
library(ggplot2)
library(lubridate)
library(tidyr)
library(sociome)
library(gtsummary)
library(flextable)
library(officer)
library(corrr)
library(GGally)
library(car)
library(CARBayesST)
library(spdep)
library(exactextractr)

setwd("C:/Users/barboza-salerno.1/Downloads/urbanization flexurba")

cbg <- block_groups(state = "IL", county = "Cook", year = 2020) %>%
  st_transform("EPSG:4326")

years <- 2020:2023
pop_list <- lapply(years, function(yr) {
  get_acs(
    geography = "block group",
    variables = "B01003_001",
    year = yr,
    survey = "acs5",
    state = "IL",
    county = "Cook",
    output = "wide"
  ) %>% 
    dplyr::select(GEOID, population = B01003_001E) %>%
    rename_with(~ paste0("pop_", yr), .cols = "population")
})

pop_df <- Reduce(function(x, y) left_join(x, y, by = "GEOID"), pop_list)

r4 <- rast("data/GHS_SMOD_E2020_GLOBE_R2023A_54009_1000_V2_0_R4_C11.tif")
r5 <- rast("data/GHS_SMOD_E2020_GLOBE_R2023A_54009_1000_V2_0_R5_C11.tif")
ghs_smod <- mosaic(r4, r5)

cbg_proj <- st_transform(cbg, crs = crs(ghs_smod))
cbg_vect <- vect(cbg_proj)

my_modal <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA)
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

urban_extract <- terra::extract(ghs_smod, cbg_vect, fun = my_modal)
cbg_proj$ghs_smod_code <- urban_extract[, 2]

ghs_labels <- c(
  "0" = "No Data", "10" = "Water", "11" = "Very Low Density Rural",
  "12" = "Low Density Rural", "13" = "Rural Cluster",
  "21" = "Suburban or Peri-Urban", "22" = "Semi-Dense Urban Cluster",
  "23" = "Dense Urban Cluster", "30" = "Urban Centre"
)
cbg_proj$urban_label <- ghs_labels[as.character(cbg_proj$ghs_smod_code)]

green_access <- read_csv("data/acs_CHI_20200622_WEDAM.csv",
                         locale = locale(decimal_mark = ".", grouping_mark = ""),
                         col_types = cols(GEOID = col_character())) %>%
  dplyr::select(GEOID, acres_c15)

overdose_df <- read_csv("data/Medical_Examiner_Case_Archive_20250614 (1).csv") %>%
  filter(!is.na(latitude), !is.na(longitude)) 

overdose_sf <- st_as_sf(overdose_df, coords = c("longitude", "latitude"), crs = 4326)
overdose_sf$year <- year(mdy_hms(overdose_sf$`Date of Death`))
overdose_sf <- overdose_sf %>% filter(year <= 2023)

overdose_joined <- st_join(overdose_sf, cbg, left = FALSE)

overdose_counts <- overdose_joined %>%
  st_drop_geometry() %>%
  group_by(GEOID, year) %>%
  summarise(overdose_n = n(), .groups = "drop") %>%
  pivot_wider(names_from = year, values_from = overdose_n, names_prefix = "overdose_") %>%
  mutate(across(starts_with("overdose_"), ~ replace_na(.x, 0)))

cbg_final <- cbg_proj %>%
  mutate(GEOID = as.character(GEOID)) %>%
  left_join(pop_df, by = "GEOID") %>%
  left_join(green_access, by = "GEOID") %>%
  left_join(overdose_counts, by = "GEOID") %>%
  mutate(across(starts_with("overdose_"), ~ replace_na(.x, 0)))

chi_adi <- get_adi(
  geography = "block group",
  state = "IL",
  county = "Cook",
  year = 2020,
  geometry = FALSE
) %>%
  left_join(cbg_final, by = "GEOID")

cbg_glm <- chi_adi %>%
  st_drop_geometry() %>%
  filter(
    !is.na(overdose_2020),
    !is.na(pop_2020),
    pop_2020 > 0
  ) %>%
  mutate(
    urban_label = factor(urban_label),
    log_pop = log(pop_2020)
  )

model_glm <- glm(
  overdose_2020 ~ acres_c15 + urban_label + ADI,
  family = poisson(link = "log"),
  offset = log_pop,
  data = cbg_glm
)

summary(model_glm)

cbg_glm_quint <- cbg_glm %>%
  filter(!is.na(acres_c15), !is.na(ADI)) %>%
  mutate(
    park_q = factor(ntile(acres_c15, 5)),
    ADI_q = factor(ntile(ADI, 5))
  )

model_quint <- glm(
  overdose_2020 ~ park_q + ADI_q + urban_label,
  family = poisson(link = "log"),
  offset = log_pop,
  data = cbg_glm_quint
)

summary(model_quint)

cbg_long <- chi_adi %>%
  st_drop_geometry() %>%
  dplyr::select(
    GEOID,
    ADI,
    urban_label,
    acres_c15,
    starts_with("overdose_"),
    starts_with("pop_")
  ) %>%
  pivot_longer(
    cols = starts_with("overdose_"),
    names_to = "year",
    names_prefix = "overdose_",
    values_to = "overdose_count"
  ) %>%
  mutate(year = as.integer(year)) %>%
  filter(
    !is.na(overdose_count),
    !is.na(pop_2020),
    pop_2020 > 0,
    !is.na(acres_c15),
    !is.na(ADI)
  ) %>%
  mutate(
    log_pop = log(pop_2020),
    ADI_q = factor(ntile(ADI, 5)),
    urban_label = factor(urban_label),
    year = factor(year)
  )

chi_adi_sf <- cbg_final %>%
  dplyr::select(GEOID, geometry) %>%
  left_join(chi_adi, by = "GEOID")

overdose_sf <- st_transform(overdose_sf, st_crs(chi_adi_sf))

library(exactextractr)

built_crop <- rast("data/built_crop_norm.tif") / 1e6
light_crop <- rast("data/light_crop.tif")
light_norm <- (light_crop - minmax(light_crop)[1]) / (minmax(light_crop)[2] - minmax(light_crop)[1])

chi_adi_sf$built_env <- exact_extract(built_crop, chi_adi_sf, "median")
chi_adi_sf$light_env <- exact_extract(light_norm, chi_adi_sf, "median")

red_band <- rast("LC09_L1TP_023031_20230722_20230723_02_T1_B4.TIF") / 10000
nir_band <- rast("LC09_L1TP_023031_20230722_20230723_02_T1_B5.TIF") / 10000

prep_to_chicago <- function(raster_input, boundary_sf) {
  boundary_sf <- st_transform(boundary_sf, crs(raster_input))
  boundary_vect <- vect(boundary_sf)
  raster_input %>% crop(boundary_vect) %>% mask(boundary_vect)
}

ndvi_raster <- prep_to_chicago((nir_band - red_band) / (nir_band + red_band), chi_adi_sf)
chi_adi_sf$ndvi_mean <- exact_extract(ndvi_raster, chi_adi_sf, "median")

final_df <- chi_adi_sf %>%
  dplyr::rename(geometry = geometry.x) %>%
  dplyr::select(
    GEOID,
    ADI,
    built_env,
    light_env,
    ndvi_mean,
    urban_label,
    acres_c15,
    starts_with("overdose_"),
    starts_with("pop_"),
    geometry
  ) %>%
  mutate(
    acres_c15_z = scale(acres_c15)[, 1],
    ndvi_mean_z = scale(ndvi_mean)[, 1]
  )

library(GGally)
library(corrr)

corr_df <- final_df %>%
  st_drop_geometry() %>%
  dplyr::select(
    ADI,
    built_env,
    light_env,
    ndvi_mean,
    acres_c15
  ) %>%
  na.omit()

cor_matrix <- cor(corr_df, use = "complete.obs")

corr_df_renamed <- corr_df %>%
  dplyr::rename(
    `Area Deprivation Index` = ADI,
    `Built-up Area Intensity` = built_env,
    `Nighttime Light Intensity` = light_env,
    `Vegetative Greenness (NDVI)` = ndvi_mean,
    `Park Access (Acres within CBG)` = acres_c15
  )

GGally::ggcorr(
  corr_df_renamed,
  label = TRUE,
  label_alpha = TRUE,
  label_round = 2,
  hjust = 0.75,
  size = 3
)

(cor_matrix <- cor(corr_df_renamed, use = "complete.obs"))

cbg_long <- final_df %>%
  st_drop_geometry() %>%
  dplyr::select(
    GEOID,
    ADI,
    built_env,
    light_env,
    ndvi_mean_z,
    acres_c15_z,
    urban_label,
    starts_with("overdose_"),
    starts_with("pop_")
  ) %>%
  pivot_longer(
    cols = starts_with("overdose_"),
    names_to = "year",
    names_prefix = "overdose_",
    values_to = "overdose_count"
  ) %>%
  mutate(
    overdose_count = replace_na(overdose_count, 0), # Note in final paper impute the mean or median except for overdose count
    pop_2020 = replace_na(pop_2020, 0),
    acres_c15_z = replace_na(acres_c15_z, 0),
    ADI = replace_na(ADI, 0),
    ndvi_mean_z = replace_na(ndvi_mean_z, 0),
    built_env = replace_na(built_env, 0),
    light_env = replace_na(light_env, 0),
    year = as.integer(year)
  ) %>%
  filter(
    pop_2020 > 0
  ) %>%
  mutate(
    log_pop = log(pop_2020),
    ADI_q = factor(ntile(ADI, 5)),
    urban_label = factor(urban_label),
    year = factor(year)
  )

model_long <- glm(
  overdose_count ~ as.numeric(year) + ADI_q + ndvi_mean_z + acres_c15_z + urban_label + built_env + light_env,
  family = poisson(link = "log"),
  offset = log_pop,
  data = cbg_long
)

summary(model_long)

library(car)

model_for_vif <- glm(
  overdose_count ~ as.numeric(year) + ADI_q + ndvi_mean_z + acres_c15_z + urban_label + built_env + light_env,
  family = poisson(link = "log"),
  data = cbg_long
)

vif(model_for_vif)

cbg_long <- cbg_long %>%
  mutate(
    urban_label_recode = dplyr::case_when(
      urban_label %in% c("Very Low Density Rural", "Low Density Rural") ~ "Rural",
      urban_label %in% c("Rural Cluster", "Semi-Dense Urban Cluster", "Suburban or Peri-Urban") ~ "Transitional",
      urban_label %in% c("Dense Urban Cluster", "Urban Centre") ~ "Urban",
      urban_label == "Water" ~ "Water",
      TRUE ~ NA_character_
    ),
    urban_label_recode = factor(urban_label_recode, levels = c("Urban", "Rural", "Transitional", "Water"))
  )

cbg_long <- cbg_long %>%
  dplyr::mutate(
    urban_label_recode = dplyr::case_when(
      urban_label %in% c("Urban Centre", "Dense Urban Cluster") ~ "Urban",
      urban_label == "Suburban or Peri-Urban" ~ "Suburban",
      urban_label == "Semi-Dense Urban Cluster" ~ "Transitional",
      urban_label %in% c("Low Density Rural", "Very Low Density Rural", "Rural Cluster") ~ "Rural",
      urban_label == "Water" ~ NA_character_
    ),
    urban_label_recode = factor(urban_label_recode, levels = c("Urban", "Suburban", "Transitional", "Rural"))
  )



model_recode <- glm(
  formula = overdose_count ~ as.numeric(year) + ADI_q + ndvi_mean_z + acres_c15_z + urban_label_recode + built_env + light_env,
  family = poisson(link = "log"),
  data = cbg_long,
  offset = log_pop
)

summary(model_recode)

library(spdep)

# Sort spatial data to match long-form data
cbg_final <- cbg_final[order(cbg_final$GEOID), ]
cbg_long <- cbg_long[order(cbg_long$GEOID, cbg_long$year), ]

cbg_final <- cbg_final %>%
  filter(GEOID %in% cbg_long$GEOID)

nb <- poly2nb(cbg_final)
W <- nb2mat(nb, style = "B")

cbg_final <- cbg_final %>%
  mutate(
    total_pop_2020 = sum(pop_2020, na.rm = TRUE),
    total_pop_2021 = sum(pop_2021, na.rm = TRUE),
    total_pop_2022 = sum(pop_2022, na.rm = TRUE),
    total_pop_2023 = sum(pop_2023, na.rm = TRUE),
    
    total_od_2020 = sum(overdose_2020, na.rm = TRUE),
    total_od_2021 = sum(overdose_2021, na.rm = TRUE),
    total_od_2022 = sum(overdose_2022, na.rm = TRUE),
    total_od_2023 = sum(overdose_2023, na.rm = TRUE),
    
    rate_2020 = total_od_2020 / total_pop_2020,
    rate_2021 = total_od_2021 / total_pop_2021,
    rate_2022 = total_od_2022 / total_pop_2022,
    rate_2023 = total_od_2023 / total_pop_2023,
    
    expected_2020 = pop_2020 * rate_2020,
    expected_2021 = pop_2021 * rate_2021,
    expected_2022 = pop_2022 * rate_2022,
    expected_2023 = pop_2023 * rate_2023
  ) %>%
  dplyr::select(-starts_with("total_"), -starts_with("rate_"))  # Cleanup

expected_long <- cbg_final %>%
  dplyr::select(GEOID, expected_2020, expected_2021, expected_2022, expected_2023) %>%
  pivot_longer(
    cols = starts_with("expected_"),
    names_to = "year",
    names_prefix = "expected_",
    values_to = "expected"
  ) %>%
  mutate(year = as.numeric(year))

cbg_long <- cbg_long %>%
  mutate(year = as.numeric(as.character(year))) %>%
  left_join(expected_long, by = c("GEOID", "year")) %>%
  mutate(
    expected = ifelse(expected == 0, 0.0001, expected),  # Avoid log(0) issues
    year_c = year - 2020  # Center year for interpretability
  )

length(unique(cbg_model$cbg_id))
names(cbg_model)
length(unique(cbg_model$GEOID))  # should be 3988
length(unique(cbg_model$year))   # already confirmed as 4

cbg_model <- cbg_model %>%
  arrange(GEOID, year)
length(cbg_model$expected) == nrow(cbg_model)

cbg_model %>%
  dplyr::count(GEOID) %>%
  dplyr::filter(n != 4)

cbg_model <- cbg_model %>%
  dplyr::arrange(GEOID, year)
library(spdep)
library(sf)

cbg_units <- cbg_model %>%
  dplyr::select(GEOID, geometry) %>%
  dplyr::distinct()

cbg_units <- sf::st_as_sf(cbg_units)
sf::st_is_valid(cbg_units) %>% table()

nb <- poly2nb(cbg_units, queen = TRUE)
W <- nb2mat(nb, style = "B", zero.policy = TRUE)
dim(W)  # Should be 3988 x 3988

library(CARBayesST)

cbg_long <- cbg_long %>%
  dplyr::filter(!is.na(urban_label_recode))

model_car <- ST.CARar(
  formula = overdose_count ~ year_c + ADI_q + ndvi_mean_z + acres_c15_z + urban_label_recode + built_env + light_env + offset(log(expected)),
  family = "poisson",
  data = cbg_long,
  W = W,
  burnin = 5000,
  n.sample = 20000,
  thin = 10,
  n.chains = 3,
  n.cores = 3,
  AR = 1,
  MALA = TRUE,
  verbose = TRUE
)

print(model_car)

fixed_effects <- model_car$summary.results

IRRs <- exp(fixed_effects[, c("Mean", "2.5%", "97.5%")])
IRRs

cbg_long$RR <- model_car$fitted.values

map_data <- cbg_long %>%
  filter(year == 2022) %>%
  left_join(cbg_final, by = "GEOID") %>%
  st_as_sf()

ggplot() +
  geom_sf(data = map_data, aes(fill = RR), color = NA) +
  scale_fill_viridis_c(name = "Relative Risk", option = "plasma") +
  theme_minimal() +
  labs(title = "Overdose Relative Risk by Census Tract, 2022")

cbg_long <- cbg_long %>%
  mutate(RR = ifelse(is.na(RR), 0, RR))

map_data <- cbg_long %>%
  left_join(cbg_final, by = "GEOID") %>%
  st_as_sf()

ggplot() +
  geom_sf(data = map_data, aes(fill = RR), color = NA) +
  scale_fill_viridis_c(name = "Relative Risk", option = "plasma", na.value = "grey90") +
  facet_wrap(~ year) +
  theme_void() +
  labs(title = "Overdose Relative Risk by Census Tract (2020–2023)")

map_data$high_risk <- map_data$RR > 1 & map_data$year == 2022

ggplot() +
  geom_sf(data = map_data %>% filter(year == 2022), aes(fill = high_risk), color = NA) +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "lightgrey"), name = "High Risk (RR > 1)") +
  theme_minimal() +
  labs(title = "High Overdose Risk Areas, 2022")

rr_trend <- cbg_long %>%
  group_by(year_c) %>%
  summarise(mean_RR = mean(RR))

ggplot(rr_trend, aes(x = year_c + 2020, y = mean_RR)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_point(color = "blue", size = 2) +
  theme_minimal() +
  labs(title = "Modeled Mean Relative Risk of Overdose by Year",
       y = "Mean Relative Risk",
       x = "Year")
obs_trend <- cbg_long %>%
  group_by(year) %>%
  summarise(total_overdoses = sum(overdose_count))

ggplot(obs_trend, aes(x = year, y = total_overdoses)) +
  geom_line(color = "darkred", linewidth = 1) +
  geom_point(color = "darkred", size = 2) +
  theme_minimal() +
  labs(title = "Observed Overdose Counts by Year",
       y = "Total Overdose Count",
       x = "Year")

library(biscale)
library(cowplot)

cookbounds <- st_read("data/cook_bounds.shp") %>%
  st_transform(26971)

chicago_bound <- st_read("data/chicago_bounds.shp")

final_df <- st_transform(final_df, 26971)

final_df <- st_intersection(final_df, cookbounds)

cbg_bi <- bi_class(
  final_df,
  x     = ADI,        # Continuous ADI
  y     = light_env,  # Nighttime light intensity
  style = "quantile",
  dim   = 4
)

map <- ggplot() +
  geom_sf(
    data  = cbg_bi,
    aes(fill = bi_class),
    color = NA,
    size  = 0.1
  ) +
  geom_sf(
    data  = chicago_bound,
    fill  = NA,
    color = "black",
    linewidth = 0.7
  ) +
  bi_scale_fill(
    pal        = "BlueOr",
    dim        = 4,
    flip_axes  = FALSE
  ) +
  theme_void() +
  labs(
    title    = "Bivariate Map: Area Deprivation Index & Nighttime Light Intensity",
    subtitle = "Higher ADI indicates greater deprivation; higher light indicates more human activity",
    caption  = "ADI (X-axis) & Nighttime Light Intensity (Y-axis)"
  ) +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    plot.caption  = element_text(size = 10, face = "bold", hjust = 1)
  )

legend <- bi_legend(
  pal    = "BlueOr",
  dim    = 4,
  xlab   = "Higher Deprivation →",
  ylab   = "Higher Light Intensity ↑",
  size   = 8
)

final_plot <- cowplot::ggdraw() +
  draw_plot(map,   0,    0, 1, 1) +
  draw_plot(legend, 0.65, 0.05, 0.3, 0.3)

final_plot

cbg_long$high_risk <- ifelse(cbg_long$RR > 1.5, 1, 0)  # Example: RR 50% above expected

cbg_long %>%
  filter(high_risk == TRUE) %>%
  group_by(year) %>%
  summarise(high_risk_count = n())

map_data <- cbg_long %>%
  dplyr::left_join(cbg_final, by = "GEOID") %>%
  st_as_sf()

ggplot() +
  geom_sf(data = map_data, aes(fill = as.factor(high_risk)), color = NA) +
  scale_fill_manual(values = c("0" = "lightgrey", "1" = "darkred"),
                    name = "High-Risk Area",
                    labels = c("Not High-Risk", "RR > 1.5")) +
  facet_wrap(~ year) +
  theme_void() +
  labs(title = "High-Risk Overdose Areas by Year",
       subtitle = "Relative Risk greater than 1, Cook County, IL",
       caption = "Derived from spatiotemporal CAR model fitted values.") +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold"))

# Calculate crude overdose rate for the county per year
total_od_2020 <- sum(final_df$overdose_2020, na.rm = TRUE)
total_pop_2020 <- sum(final_df$pop_2020, na.rm = TRUE)
county_rate_2020 <- total_od_2020 / total_pop_2020

total_od_2021 <- sum(final_df$overdose_2021, na.rm = TRUE)
total_pop_2021 <- sum(final_df$pop_2021, na.rm = TRUE)
county_rate_2021 <- total_od_2021 / total_pop_2021

total_od_2022 <- sum(final_df$overdose_2022, na.rm = TRUE)
total_pop_2022 <- sum(final_df$pop_2022, na.rm = TRUE)
county_rate_2022 <- total_od_2022 / total_pop_2022

total_od_2023 <- sum(final_df$overdose_2023, na.rm = TRUE)
total_pop_2023 <- sum(final_df$pop_2023, na.rm = TRUE)
county_rate_2023 <- total_od_2023 / total_pop_2023

# Calculate expected counts for each CBG
final_df <- final_df %>%
  mutate(
    expected_2020 = pop_2020 * county_rate_2020,
    expected_2021 = pop_2021 * county_rate_2021,
    expected_2022 = pop_2022 * county_rate_2022,
    expected_2023 = pop_2023 * county_rate_2023,
    
    SIR_2020 = overdose_2020 / expected_2020,
    SIR_2021 = overdose_2021 / expected_2021,
    SIR_2022 = overdose_2022 / expected_2022,
    SIR_2023 = overdose_2023 / expected_2023
  )

final_df <- final_df %>%
  mutate(
    expected_2020 = ifelse(expected_2020 == 0, 0.0001, expected_2020),
    expected_2021 = ifelse(expected_2021 == 0, 0.0001, expected_2021),
    expected_2022 = ifelse(expected_2022 == 0, 0.0001, expected_2022),
    expected_2023 = ifelse(expected_2023 == 0, 0.0001, expected_2023),
    
    SIR_2020 = overdose_2020 / expected_2020,
    SIR_2021 = overdose_2021 / expected_2021,
    SIR_2022 = overdose_2022 / expected_2022,
    SIR_2023 = overdose_2023 / expected_2023
  )

# Join decedent data to CBG-level attributes
decedent_data <- overdose_joined %>% st_drop_geometry() %>%
  left_join(
    final_df %>% dplyr::select(GEOID, ADI, built_env, light_env, ndvi_mean, acres_c15, SIR_2020, SIR_2021, SIR_2022, SIR_2023),
    by = "GEOID"
  )

decedent_data <- decedent_data %>%
  mutate(
    race_ethnicity = case_when(
      Latino == TRUE ~ "Latino",
      Race == "White" & Latino == FALSE ~ "Non-Hispanic White",
      Race == "Black" ~ "Non-Hispanic Black",
      Race == "Asian" ~ "Non-Hispanic Asian",
      Race == "Am. Indian" ~ "Non-Hispanic Indigenous",
      Race == "Other" ~ "Other",
      Race == "Unknown" ~ "Unknown",
      TRUE ~ "Missing"
    )
  )
decedent_summary <- decedent_data %>%
  st_drop_geometry() %>%
  dplyr::select(GEOID, age = Age, gender = Gender, race_ethnicity, ADI:overdose_2023) %>%
  summarise(
    mean_age = mean(age, na.rm = TRUE),
    n_male = sum(gender == "Male", na.rm = TRUE),
    n_female = sum(gender == "Female", na.rm = TRUE),
    n_other_gender = sum(!gender %in% c("Male", "Female") & !is.na(gender)),
    
    n_latino = sum(race_ethnicity == "Latino", na.rm = TRUE),
    n_nh_white = sum(race_ethnicity == "Non-Hispanic White", na.rm = TRUE),
    n_nh_black = sum(race_ethnicity == "Non-Hispanic Black", na.rm = TRUE),
    n_nh_asian = sum(race_ethnicity == "Non-Hispanic Asian", na.rm = TRUE),
    n_indigenous = sum(race_ethnicity == "Non-Hispanic Indigenous", na.rm = TRUE),
    n_other_race = sum(race_ethnicity == "Other", na.rm = TRUE),
    n_unknown_race = sum(race_ethnicity == "Unknown", na.rm = TRUE),
    
    across(ADI:overdose_2023, ~mean(.x, na.rm = TRUE), .names = "{.col}_mean")
  )
##############################
fitted_vals <- model_car$fitted.values  # matrix: rows = observations, cols = chains or summaries
cbg_long$fitted <- model_car$fitted.values
library(dplyr)
library(broom)

# Linear trend
slope_trends <- cbg_long %>%
  group_by(GEOID) %>%
  do({
    fit <- lm(fitted ~ year_c, data = .)
    tidy(fit) %>% filter(term == "year_c")
  }) %>%
  rename(slope = estimate, p_slope = p.value)

# Quadratic trend
curve_trends <- cbg_long %>%
  group_by(GEOID) %>%
  do({
    fit <- lm(fitted ~ year_c + I(year_c^2), data = .)
    tidy(fit) %>% filter(term == "I(year_c^2)")
  }) %>%
  rename(curvature = estimate, p_curve = p.value)

# Combine and classify
trend_summary <- slope_trends %>%
  left_join(curve_trends, by = "GEOID") %>%
  mutate(
    trend_type = case_when(
      slope > 0 & p_slope < 0.05 & curvature > 0 & p_curve < 0.05 ~ "sharply increasing",
      slope > 0 & p_slope < 0.05 ~ "increasing",
      slope < 0 & p_slope < 0.05 ~ "decreasing",
      TRUE ~ "no change"
    )
  )

cbg_classified <- cbg_long %>%
  left_join(trend_summary %>% dplyr::select(GEOID, trend_type), by = "GEOID")

# Extract distinct GEOID and urban label from the full data
urban_lookup <- cbg_long %>%
  dplyr::select(GEOID, urban_label_recode) %>%
  distinct()

# Join to trend summary
trend_summary_labeled <- trend_summary %>%
  left_join(urban_lookup, by = "GEOID")

urban_trend_summary <- trend_summary_labeled %>%
  group_by(urban_label_recode, trend_type) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(urban_label_recode) %>%
  mutate(
    total = sum(n),
    proportion = round(n / total, 3)
  ) %>%
  arrange(urban_label_recode, desc(proportion))

ggplot(urban_trend_summary, aes(x = urban_label_recode, y = proportion, fill = trend_type)) +
  geom_col(position = "stack") +
  labs(
    title = "Overdose Risk Trend Classification by Urban Category",
    x = "Urban Classification",
    y = "Proportion of GEOIDs",
    fill = "Trend Type"
  ) +
  theme_minimal()
