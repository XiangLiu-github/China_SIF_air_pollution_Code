library(tidyverse)
library(readxl)
library(tidymodels)
library(dtplyr)
library(qs)
library(fixest)
library(collapse)
library(furrr)
# source("script/loadPackages.R")
source("script/loadFunctions.R")
source("script/loadFormulas.R")

# f1 ----------------------------------------------------------------------

data <- get_data()

data <- data %>%
  group_by(crop, x, y) %>%
  mutate(vpd_mean = fmean(vpd_mean)) %>%
  group_by(crop_parent) %>%
  mutate(vpd_mean = DescTools::Winsorize(vpd_mean, minval = 3, maxval = 15)) %>%
  ungroup()

n <- 500

data <- data %>% nest(fdata = -crop_parent)

data$coefs <- map2(data$crop_parent, data$fdata, function(acrop, adata) {
  adata <- adata %>%
    mutate(AOD2 = AOD^2)

  adata <- map(c(
    "AOD", "AOD2"
  ), function(avar) {
    i(adata$region, adata %>% pull(!!avar)) %>%
      as_tibble() %>%
      `names<-`(str_c(names(.), "_", avar))
  }) %>%
    bind_cols() %>%
    bind_cols(adata)

  adata <- adata %>%
    nest(.by = county, .key = "ffdata")

  plan(multisession, workers = 4)

  res <- future_map_dfr(1:n, function(anum) {
    set.seed(anum)
    aadata <- adata %>%
      slice_sample(prop = 1, replace = T) %>%
      unnest(ffdata)

    if (acrop == "Rice") {
      model <- feols(fml_inter_full_rice, aadata,
        weights = ~fraction,
        nthreads = 6, notes = F, lean = T
      )
    } else {
      model <- feols(fml_inter_full, aadata,
        weights = ~fraction,
        nthreads = 6, notes = F, lean = T
      )
    }

    broom::tidy(model)
  }, .id = "id", .progress = T, .options = furrr_options(seed = T))

  plan(sequential)

  return(res)
})

saveRDS(data %>% select(-fdata), "data/boots_f1_hete_full.rds")
