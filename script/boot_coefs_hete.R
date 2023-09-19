library(tidyverse)
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

data <- data %>%
  mutate(coefs = map2(crop_parent, fdata, function(acrop, adata) {
    adata <- adata %>%
      nest(.by = county, .key = "ffdata")

    plan(multisession, workers = 4)

    res <- future_map_dfr(1:n, function(anum) {
      set.seed(anum)
      aadata <- adata %>%
        slice_sample(prop = 1, replace = T) %>%
        unnest(ffdata)

      if (acrop == "Rice") {
        model <- feols(fml_inter_v, aadata,
          weights = ~fraction,
          nthreads = 6, notes = F, lean = T
        )
      } else {
        model <- feols(fml_inter_all, aadata,
          weights = ~fraction,
          nthreads = 6, notes = F, lean = T
        )
      }

      broom::tidy(model)
    }, .id = "id", .progress = T, .options = furrr_options(seed = T))

    plan(sequential)

    return(res)
  }))

saveRDS(data %>% select(-fdata), "data/boots_f1_hete.rds")
