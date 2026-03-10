##################################################################
#                               
# ESM 211: Applied Population Ecology
# 
# Project: Effects of Wind Installation on Golden Eagle Populations
# 
# Version: 1
# Feb 20, 2026           
# Ben Ladabaum and Scott Schwartz                 
#                               
# Description:  Build dataset and analyze effects of wind installations 
#               on Golden Eagles in CA
#
# Notes:    
# 
#
##################################################################

# Clean the environment
rm(list=ls())

# Load required packages
library(here)
library(janitor)
library(tidyverse)
library(ggplot2)
library(sf)
library(writexl)
library(readxl)
library(lme4) # for negative binomial and poisson regressions
library(tigris) # for adding CA to map
library(broom.mixed) # for outputting regression
library(knitr) # for outputting regression

# load data
routes = read_csv(here("Data", "routes.csv")) |>
  clean_names() |>
  filter(state_num=="14") |> # filter to CA
  filter(route_type_id==1) 

ca_golden_eagles = read_csv(here("Data", "Califor.csv")) |>
  clean_names() |>
  filter(aou=="03490") |> # filter to golden eagles
  left_join(y = routes, by = c("route"))

wind_data = read_csv(here("Data", "uswtdb_V8_2_20251210.csv")) |>
  filter(t_state=="CA") |> # filter to CA 
  group_by(p_name, p_year, p_tnum, p_cap) |> # calculate mean turbine location for each project
  summarise(
    latitude  = mean(ylat,  na.rm = TRUE),
    longitude = mean(xlong, na.rm = TRUE)
  ) |>
  filter(!is.na(p_year)) # filter out obs with unknown year

ge_survey = read_excel(here("Data", "midwinter_bald_eagle_survey_appendix1.xlsx")) |>
  clean_names() |>
  filter(state=="CA")


# prepare data for spacial merge
wind_sf <- wind_data |>
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) |>
  st_transform(crs = 5070) # Transform to a projected CRS for distance-based operations (e.g., meters)

golden_eagles_sf <- ca_golden_eagles |>
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) |>
  st_transform(crs = 5070)

ge_survey_sf = ge_survey |>
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) |>
  st_transform(crs = 5070)


# Spatial join — find bird routes within buffer of each wind project centroid
######
# Good place for sensitivity analyses (change buffer distance and see how this matters)
######
#wind_buffered <- st_buffer(wind_sf, dist = 50000)  # 50,000 meter radius around each project
wind_buffered <- st_buffer(wind_sf, dist = 100000) # 100,000 meter radius around each project

golden_eagles_w_wind = st_join(golden_eagles_sf, 
                               wind_buffered, join = st_intersects, left = TRUE)

ge_survey_w_wind = st_join(ge_survey_sf, 
                        wind_buffered, join = st_intersects, left = TRUE)

# export data (not necessary for full program; just to be able to look at data better)
golden_eagle_routes_w_wind =
  st_drop_geometry(golden_eagles_w_wind) 
  
write_xlsx(golden_eagle_routes_w_wind, here("Output", "golden_eagle_routes_w_wind.xlsx")) 


# Define pre and post treatment
golden_eagles_w_wind = golden_eagles_w_wind |>
  mutate(treatment = case_when(
    is.na(p_name) ~ 0,
    !is.na(p_name) ~1
  )) |>
  mutate(post = case_when(
    treatment==0 ~ 0,
    treatment==1 & year<p_year ~ 0,
    treatment==1 & year>=p_year ~ 1
  )) |>
  mutate(treatment = as.factor(treatment)) |>
  mutate(year_c = year - min(year)) #baseline year as years since first obs; important so neg binom regression doesn't blow up

ge_survey_w_wind = ge_survey_w_wind |>
  mutate(treatment = case_when(
    is.na(p_name) ~ 0,
    !is.na(p_name) ~1
  )) |>
  mutate(post = case_when(
    treatment==0 ~ 0,
    treatment==1 & year<p_year ~ 0,
    treatment==1 & year>=p_year ~ 1
  )) |>
  mutate(treatment = as.factor(treatment)) |>
  mutate(year_c = year - min(year)) #baseline year as years since first obs; important so neg binom regression doesn't blow up



### Maps
# create layer for CA
california <- states(cb = TRUE) |>
  filter(NAME == "California") |>
  st_transform(st_crs(golden_eagles_sf))  # match CRS to your existing data

ggplot() +
  geom_sf(data = california, fill = "grey95", color = "black", linewidth = 0.5) +
  geom_sf(data = wind_buffered, aes(fill = "100km Buffer"), color = "blue",
          alpha = 0.3, linewidth = 0.5) +
  geom_sf(data = golden_eagles_sf, aes(color = "Eagle Routes"),
          alpha = 0.5, linewidth = 0.5) +
  geom_sf(data = wind_sf, aes(color = "Wind Projects"),
          size = 2, shape = 17) +
  scale_color_manual(name = NULL,
                     values = c("Eagle Routes" = "orange", 
                                "Wind Projects" = "blue")) +
  scale_fill_manual(name = NULL,
                    values = c("100km Buffer" = "lightblue")) +
  labs(title = "Golden Eagle Routes and Wind Project Locations: Bird Breeding Survey") +
  #coord_sf(datum = NA) +
  theme_minimal() +
  theme(legend.position = "bottom")  # move legend wherever you like: "top", "left", "right"

ggplot() +
  geom_sf(data = california, fill = "grey95", color = "black", linewidth = 0.5) +
  geom_sf(data = wind_buffered, aes(fill = "100km Buffer"), color = "blue",
          alpha = 0.3, linewidth = 0.5) +
  geom_sf(data = ge_survey_sf, aes(color = "Eagle Routes"),
          alpha = 0.5, linewidth = 0.5) +
  geom_sf(data = wind_sf, aes(color = "Wind Projects"),
          size = 2, shape = 17) +
  scale_color_manual(name = NULL,
                     values = c("Eagle Routes" = "red", 
                                "Wind Projects" = "blue")) +
  scale_fill_manual(name = NULL,
                    values = c("100km Buffer" = "lightblue")) +
  labs(title = "Golden Eagle Routes and Wind Project Locations: Golden Eagle Survey") +
  #coord_sf(datum = NA) +
  theme_minimal() +
  theme(legend.position = "bottom")  # move legend wherever you like: "top", "left", "right"



### run models

# Consistency checks - greater mean than variance could mean under-dispersion
mean(golden_eagles_w_wind$species_total) #~1.51
var(golden_eagles_w_wind$species_total) #0.95

mean(ge_survey_w_wind$totalge) #1.04
var(ge_survey_w_wind$totalge) #14.35

# 1) Using "post" as key variable - estimates avg treatment effect (discrete before/after effect)
neg_binom_model <- glmer.nb(
  species_total ~ year_c + treatment
                + post               # key coefficient
                + (1 | route),       # random effects for route
  data = golden_eagles_w_wind
)
summary(neg_binom_model)


neg_binom_model_gesurvey <- glmer.nb(
  totalge ~ year_c + treatment
  + post                     # key coefficient
  + (1 | sitenumb),          # random intercept for site
  data = ge_survey_w_wind
)
summary(neg_binom_model_gesurvey)




# poisson - if no overdispersion and mean = variance
model_poisson = glmer(
  species_total ~ year_c + treatment + post + (1 | route),
  data   = golden_eagles_w_wind,
  family = poisson
)
summary(model_poisson)

model_poisson_gesurvey = glmer(
  totalge ~ year_c + treatment + post + (1 | sitenumb),
  data   = ge_survey_w_wind,
  family = poisson
)
summary(model_poisson_gesurvey)



## Output regression results
tidy(model_poisson, effects = "fixed", conf.int = TRUE) |>
  kable(digits = 3)

tidy(neg_binom_model_gesurvey, effects = "fixed", conf.int = TRUE) |>
  kable(digits = 3)





# Model assessment for bird breeding survey
AIC(neg_binom_model, model_poisson) # better for predictions
BIC(neg_binom_model, model_poisson) # better for causal anlaysis

# Model assessment eagle survey
AIC(neg_binom_model_gesurvey, model_poisson_gesurvey) # better for predictions
BIC(neg_binom_model_gesurvey, model_poisson_gesurvey) # better for causal anlaysis






