library(terra)
library(data.table)
library(tidyverse)
library(tidymodels)
library(dtplyr)
library(qs)
library(matrixStats)
library(furrr)
library(collapse)
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

mvpd <- read_rds("data/outputs/vpd/tidied.rds") %>%
  lazy_dt() %>%
  group_by(x, y, month) %>%
  summarise(vpd_mean = mean(vpd)) %>%
  ungroup() %>%
  as_tibble()

mvpd <- cldr %>%
  lazy_dt() %>%
  inner_join(mvpd) %>%
  group_by(crop, x, y, year) %>%
  summarise(vpd_mean = mean(vpd_mean), .groups = "drop") %>%
  as_tibble() %>%
  group_by(crop, x, y) %>%
  mutate(
    vpd_mean = fmean(vpd_mean),
    crop_parent = fifelse(str_detect(crop, "Rice"), "Rice", crop)
  ) %>%
  group_by(crop_parent) %>%
  mutate(vpd_mean = DescTools::Winsorize(vpd_mean, probs = c(0.01, 0.99))) %>%
  ungroup() %>%
  select(-crop_parent) %>%
  nest(vdata = -c(crop, year)) %>%
  arrange(crop, year)

f1 <- read_rds("data/boots_f1_hete_vpd.rds") %>%
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

plan(multisession, workers = 4)

impacts <-
  reduce(list(AOD, mvpd, f1, aod_ctr), inner_join)

# acrop = 'Maize'
# adata = impacts$fdata[[1]]
# avdata = impacts$vdata[[1]]
# coefs = impacts$coefs[[1]]
# aod_data = impacts$aod_data[[1]]

impacts <- impacts %>%
  mutate(AOD_results = future_pmap(list(crop, fdata, vdata, coefs, aod_data), function(acrop, adata, avdata, coefs, aod_data) {
    coefs <- coefs %>%
      select(id, term, estimate) %>%
      pivot_wider(names_from = term, values_from = estimate) %>%
      select(contains("AOD"))

    cal_data <- adata %>%
      inner_join(aod_data, by = join_by(x, y)) %>%
      inner_join(avdata, by = join_by(x, y)) %>%
      mutate(AOD_cft = fifelse(AOD_cft > AOD, AOD, AOD_cft))

    rel_X <- cal_data %>%
      mutate(
        AOD1 = AOD_cft - AOD, AOD2 = AOD_cft^2 - AOD^2,
        vAOD1 = AOD1 * vpd_mean, vAOD2 = AOD2 * vpd_mean
      )

    rel_X <- rel_X %>% select(AOD1, AOD2, vAOD1, vAOD2)

    stopifnot(sum(rel_X$AOD1 > 0) == 0)

    # cal relative impact
    rel_results <- tcrossprod(as.matrix(rel_X), as.matrix(coefs)) %>%
      expm1() %>%
      rowQuantiles(probs = c(0.05, 0.5, 0.95)) %>%
      as_tibble() %>%
      bind_cols(cal_data %>% select(x, y, fraction, pm_level))

    return(rel_results)
  }, .progress = T)) %>%
  select(-c(fdata, vdata, coefs, aod_data))

plan(sequential)

qsave(impacts, "data/impacts_aerosol_hete_vpd.qs", nthread = qn)

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
  summarise(AOT40 = sum(AOT40), .groups = "drop") %>%
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
      summarise(AOT40_cft = sum(AOT40_cft), .groups = "drop") %>%
      as_tibble()
  }, .progress = T)) %>%
  unnest() %>%
  nest(ozone_data = -c(crop, year))

impacts <-
  reduce(list(ozone, f1, mvpd, ozone_ctr), inner_join)

# acrop = 'Rice(SR&ER)'
# adata = impacts$fdata[[1]]
# avdata = impacts$vdata[[1]]
# coefs = impacts$coefs[[1]]
# ozone_data = impacts$ozone_data[[1]]

plan(multisession, workers = 4)

impacts <- impacts %>%
  mutate(ozone_results = future_pmap(list(crop, fdata, vdata, coefs, ozone_data), function(acrop, adata, avdata, coefs, ozone_data) {
    coefs <- coefs %>%
      select(id, term, estimate) %>%
      pivot_wider(names_from = term, values_from = estimate) %>%
      select(contains("AOT40"))

    cal_data <- adata %>%
      inner_join(ozone_data, by = c("x", "y")) %>%
      inner_join(avdata, by = join_by(x, y)) %>%
      mutate(AOT40_cft = fifelse(AOT40_cft > AOT40, AOT40, AOT40_cft))

    rel_X <- cal_data %>%
      mutate(
        AOT401 = AOT40_cft - AOT40,
        vAOT401 = AOT401 * vpd_mean
      )

    rel_X <- rel_X %>% select(AOT401, vAOT401)

    stopifnot(sum(rel_X$AOT401 > 0) == 0)

    # cal relative impact
    rel_results <- tcrossprod(as.matrix(rel_X), as.matrix(coefs)) %>%
      expm1() %>%
      rowQuantiles(probs = c(0.05, 0.5, 0.95)) %>%
      as_tibble() %>%
      bind_cols(cal_data %>% select(x, y, fraction, peak_level))

    return(rel_results)
  }, .progress = T)) %>%
  select(-c(fdata, vdata, coefs, ozone_data))

plan(sequential)

qsave(impacts, "data/impacts_ozone_AOT40_hete_vpd.qs", nthread = qn)

# impacts$ozone_results[[60]] %>%
#   filter(peak_level == 30) %>%
#   select(x, y, `50%`) %>%
#   rast(type = "xyz") %>%
#   plot()
