# source("script/loadPackages.R")
library(terra)
library(data.table)
library(tidyverse)
library(tidymodels)
library(dtplyr)
library(qs)
library(matrixStats)
library(furrr)
source("script/loadFunctions.R")

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
        pivot_longer(
          -c(x, y),
          names_to = "year",
          names_transform = list(year = as.integer),
          values_drop_na = T,
          values_to = "fraction"
        )
    })
  ) %>%
  unnest() %>%
  trim_xy()

f1 <- read_rds("data/boots_f1_w126.rds") %>%
  mutate(crop = crop_parent) %>%
  bind_rows((.) %>% filter(crop == "Rice") %>% mutate(crop = "Rice(LR)")) %>%
  mutate(crop = fifelse(crop == "Rice", "Rice(SR&ER)", crop))

# aerosol -----------------------------------------------------------------

AOD <- read_rds("data/outputs/aerosol/AOD.rds")

AOD <- cldr %>%
  lazy_dt() %>%
  inner_join(AOD) %>%
  group_by(crop, x, y, year) %>%
  summarise(AOD = mean(AOD), .groups = "drop") %>%
  inner_join(fraction) %>%
  as_tibble() %>%
  nest(fdata = -c(crop, year)) %>%
  arrange(crop, year)

aod_ctr <- read_rds("data/aod_cft.rds")

aod_ctr <- aod_ctr %>%
  mutate(aod_data = map(aod_data, function(adata) {
    cldr %>%
      lazy_dt() %>%
      inner_join(adata, by = c("x", "y", "month")) %>%
      group_by(crop, year, x, y) %>%
      summarise(AOD_cft = mean(AOD_cft), .groups = "drop") %>%
      as_tibble()
  }, .progress = T)) %>%
  unnest() %>%
  nest(aod_data = -c(crop, year))

impacts <-
  reduce(list(AOD, f1, aod_ctr), inner_join)

# adata = impacts$fdata[[71]]
# coefs = impacts$coefs[[71]]
# aod_data = impacts$aod_data[[71]]

plan(multisession, workers = 6)

impacts <- impacts %>%
  mutate(AOD_results = future_pmap(list(fdata, coefs, aod_data), function(adata, coefs, aod_data) {
    coefs <- coefs %>%
      select(id, term, estimate) %>%
      pivot_wider(names_from = term, values_from = estimate) %>%
      select(contains("AOD"))

    cal_data <- adata %>%
      inner_join(aod_data, by = c("x", "y")) %>%
      mutate(AOD_cft = fifelse(AOD_cft > AOD, AOD, AOD_cft))

    rel_X <- cal_data %>%
      mutate(AOD1 = AOD_cft - AOD, AOD2 = AOD_cft^2 - AOD^2) %>%
      select(AOD1, AOD2)

    stopifnot(sum(rel_X$AOD1 > 0) == 0)

    # cal relative impact
    rel_results <- tcrossprod(as.matrix(rel_X), as.matrix(coefs)) %>%
      expm1() %>%
      rowQuantiles(probs = c(0.05, 0.5, 0.95)) %>%
      as_tibble() %>%
      bind_cols(cal_data %>% select(x, y, fraction, pm_level))

    return(rel_results)
  }, .progress = T)) %>%
  select(-c(fdata, coefs, aod_data))

plan(sequential)

qsave(impacts, "data/impacts_aerosol_W126.qs", nthread = qn)

# impacts$AOD_results[[20]] %>%
#   filter(pm_level == 50) %>%
#   select(x, y, `50%`) %>%
#   rast(type = "xyz") %>%
#   plot()

# ozone -------------------------------------------------------------------

cldr <-
  bind_rows(
    cldr %>%
      filter(crop == "Maize" & month - `GR&EM` != 6),
    cldr %>%
      filter(str_detect(crop, "Rice")),
    cldr %>%
      filter(crop == "Wheat" & MA != month)
  )

ozone <- read_rds("data/outputs/ozone/tidied.rds")

ozone <- cldr %>%
  lazy_dt() %>%
  inner_join(ozone) %>%
  group_by(crop, x, y, year) %>%
  summarise(W126 = sum(W126), .groups = "drop") %>%
  inner_join(fraction) %>%
  as_tibble() %>%
  nest(fdata = -c(crop, year)) %>%
  arrange(crop, year)

ozone_ctr <- read_rds("data/ozone_cft.rds")

ozone_ctr <- ozone_ctr %>%
  mutate(ozone_data = map(ozone_data, function(adata) {
    cldr %>%
      lazy_dt() %>%
      inner_join(adata, by = c("x", "y", "month")) %>%
      group_by(crop, year, x, y) %>%
      summarise(W126_cft = sum(W126_cft), .groups = "drop") %>%
      as_tibble()
  }, .progress = T)) %>%
  unnest() %>%
  nest(ozone_data = -c(crop, year))

impacts <-
  reduce(list(ozone, f1, ozone_ctr), inner_join)

# adata = impacts$fdata[[1]]
# coefs = impacts$coefs[[1]]
# ozone_data = impacts$ozone_data[[1]]

plan(multisession, workers = 4)

impacts <- impacts %>%
  mutate(ozone_results = future_pmap(list(fdata, coefs, ozone_data), function(adata, coefs, ozone_data) {
    coefs <- coefs %>%
      select(id, term, estimate) %>%
      pivot_wider(names_from = term, values_from = estimate) %>%
      select(contains("W126"))

    cal_data <- adata %>%
      inner_join(ozone_data, by = c("x", "y")) %>%
      mutate(W126_cft = fifelse(W126_cft > W126, W126, W126_cft))

    rel_X <- cal_data %>%
      mutate(W1261 = W126_cft - W126) %>%
      select(W1261)

    stopifnot(sum(rel_X$W1261 > 0) == 0)

    # cal relative impact
    rel_results <- tcrossprod(as.matrix(rel_X), as.matrix(coefs)) %>%
      expm1() %>%
      rowQuantiles(probs = c(0.05, 0.5, 0.95)) %>%
      as_tibble() %>%
      bind_cols(cal_data %>% select(x, y, fraction, peak_level))

    return(rel_results)
  }, .progress = T)) %>%
  select(-c(fdata, coefs, ozone_data))

plan(sequential)

qsave(impacts, "data/impacts_ozone_W126.qs", nthread = qn)

# impacts$ozone_results[[60]] %>%
#   filter(peak_level == 30) %>%
#   select(x, y, `50%`) %>%
#   rast(type = "xyz") %>%
#   plot()
