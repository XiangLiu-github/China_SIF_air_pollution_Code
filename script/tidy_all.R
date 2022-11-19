#### in FLEX
library(collapse, exclude = "F")
library(tidyverse)
library(terra)
library(Matrix)
library(dtplyr)
library(sf)
sf_use_s2(F)
library(qs)

trim_xy <- function(adata_xy) {
  adata_xy %>%
    mutate(across(c(x, y), ~ round(.x, 4)))
}

cldr <- read_rds("data/outputs/calendar/tidied.rds") %>%
  filter((MA - `GR&EM`) >= 2) %>%
  trim_xy()


fraction <-
  tibble(
    crop = c("Maize", "Rice(LR)", "Rice(SR&ER)", "Wheat"),
    data = map(c("Maize", "Rice", "Rice", "Wheat"), function(acrop) {
      afile <- str_c("data/outputs/masks/mask_", acrop, ".tif")
      rast(afile) %>%
        as.data.frame(xy = T, na.rm = F) %>%
        pivot_longer(-c(x, y),
          names_to = "year", names_transform = list(year = as.integer),
          values_drop_na = T, values_to = "fraction"
        )
    })
  ) %>%
  unnest() %>%
  trim_xy()

temp_bins <- qread("data/outputs/temp/temp_bins.qs", nthreads = 10) %>%
  unnest() %>%
  rename_with(~ str_replace(.x, "-", "_")) %>% # 2000-2021
  mutate(
    bin0 = sum(
      bin_10, bin_9, bin_8, bin_7, bin_6,
      bin_5, bin_4, bin_3, bin_2, bin_1, bin0
    ),
    bin39 = sum(
      bin39, bin40, bin41, bin42, bin43, bin44,
      bin45, bin46, bin47, bin48, bin49, bin50
    ),
    .keep = "unused"
  ) %>%
  mutate(across(contains("bin"), ~ . / 24))

df <- 4

data_bin <- temp_bins %>%
  select(num_range("bin", 0:39))

bins <- seq(0, 39, by = 1)

steps <- seq(min(bins), max(bins), df)
idx <- rep(1:length(steps), each = df)
idx <- idx[1:length(bins)]
B <- sparseMatrix(i = 1:length(bins), j = idx, x = 1)
bindata <- data_bin[, grepl("bin", names(data_bin))]
names(bindata) <- paste0("bin", bins)
bdata <- as.matrix(bindata) %*% B
bdata <- as.matrix(bdata) %>% as.data.frame()
colnames(bdata) <- paste0("step", unique(idx))

temp_bins <- temp_bins %>% bind_cols(bdata)

GOSIF <- read_rds("data/outputs/SIF/GOSIF.rds") %>%
  mutate(GOSIF = GOSIF * 0.0001) # 2000-2021
WenSIF <- read_rds("data/outputs/SIF/Wenetal.rds") # 2002-2018
RTSIF <- read_rds("data/outputs/SIF/RTSIF.rds") # 2001-2020
CSIF <- read_rds("data/outputs/SIF/CSIF.rds") # 2000-2020
tmax <- read_rds("data/outputs/temp/tmax.rds") %>%
  mutate(maxtmp = maxtmp - 273.15) # 2000-2019
prep <- read_rds("data/outputs/prep/tidied.rds") # 2000-2019
surface <- read_rds("data/outputs/soil_moisture/surface.rds") # 2000-2020
root <- read_rds("data/outputs/soil_moisture/root.rds") # 2000-2020
sm <- read_rds("data/outputs/soil_moisture/mean.rds") # 2000-2020
cloud <- read_rds("data/outputs/cloud/tidied.rds") # 2002-2020
AOD <- read_rds("data/outputs/aerosol/AOD.rds") # 2000-2020
PM25 <- read_rds("data/outputs/aerosol/PM25.rds") # 2000-2020
PM10 <- read_rds("data/outputs/aerosol/PM10.rds") # 2000-2020
ozone <- read_rds("data/outputs/ozone/tidied.rds") # 2005-2019
irrigation <- read_rds("data/outputs/irrigation/tidied.rds") # %>% # 2000-2019
# bind_rows((.) %>% filter(year == 2019) %>% mutate(year = 2020))
radiation <- read_rds("data/outputs/radiation/tidied.rds") # 2007-2018

data <-
  cldr %>%
  inner_join(fraction) %>%
  left_join(GOSIF) %>%
  left_join(WenSIF) %>%
  left_join(RTSIF) %>%
  left_join(CSIF) %>%
  left_join(temp_bins) %>%
  left_join(tmax) %>%
  left_join(prep) %>%
  left_join(surface) %>%
  left_join(root) %>%
  left_join(sm) %>%
  left_join(cloud) %>%
  left_join(AOD) %>%
  left_join(PM25) %>%
  left_join(PM10) %>%
  left_join(radiation)

shp_file <-
  st_read("data/inputs/shp/2020年县.shp") %>%
  transmute(
    county = 县级, county_code = 县级码,
    city = 地级, city_code = 地级码,
    province = 省级, province_code = 省级码
  )

xy_info <-
  data %>%
  distinct(x, y) %>%
  st_as_sf(coords = c("x", "y"), crs = st_crs(shp_file), remove = F) %>%
  st_intersection(shp_file) %>%
  st_drop_geometry()

ozone_s <-
  map(c(1:7), function(anum) {
    cldr %>%
      filter((MA - month) < anum) %>%
      inner_join(fraction) %>%
      left_join(ozone) %>%
      lazy_dt() %>%
      group_by(x, y, year, crop) %>%
      summarise(O3 = mean(O3), across(c(W126, AOT40), sum)) %>%
      as_tibble() %>%
      rename_with(~ str_c(.x, "_", anum), c(O3, W126, AOT40))
  }) %>%
  reduce(full_join)

# diff_sif = map(list(GOSIF, RTSIF, CSIF),
#                function(adata){
#                  adata %>%
#                    fgroup_by(year, x, y) %>%
#                    fdiff(t = month) %>%
#                    fungroup() %>%
#                    right_join(cldr) %>%
#                    filter(if_any(ends_with('SIF'), ~ .x >= 0)) %>%
#                    group_by(crop, x, y, year) %>%
#                    summarise(across(ends_with('SIF'), sum, .names = '{.col}_diff'), .groups = 'drop')
#                }) %>%
#   reduce(full_join)

fnl_data <-
  data %>%
  lazy_dt() %>%
  group_by(x, y, year, crop) %>%
  summarise(
    across(c(
      GOSIF, Wenetal, RTSIF, CSIF
    ), max, .names = "{.col}_peak"),
    across(c(
      GOSIF, Wenetal, RTSIF, CSIF
    ), sum, .names = "{.col}_sum"),
    across(c(
      maxtmp, cloud, AOD, PM25, PM10, 
      HE, MA, `GR&EM`,
      GOSIF, Wenetal, RTSIF, CSIF,
      fraction, DR, GR, surface, root
    ), mean),
    across(c(
      prep,
      starts_with("bin"), starts_with("step"),
      starts_with("surface_"), starts_with("root_")
    ), sum)
  ) %>%
  mutate(x_y = str_c(x, "_", y)) %>%
  inner_join(xy_info) %>%
  inner_join(irrigation) %>%
  left_join(ozone_s) %>%
  # left_join(diff_sif) %>%
  as_tibble()

qsave(fnl_data, "data/tidied.qs", preset = "fast", nthreads = 10)
