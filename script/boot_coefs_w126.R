library(tidyverse)
library(tidymodels)
library(dtplyr)
library(qs)
library(fixest)
# source("script/loadPackages.R")
source("script/loadFunctions.R")
source("script/loadFormulas.R")

# f1 ----------------------------------------------------------------------

data <- get_data()
n <- 500

data <- data %>% nest(fdata = -crop_parent)

data <- data %>%
  mutate(coefs = map(fdata, function(adata) {
    adata <- adata %>%
      nest(.by = county, .key = "ffdata")

    map_dfr(1:n, function(anum) {
      set.seed(anum)
      aadata <- adata %>%
        slice_sample(prop = 1, replace = T) %>%
        unnest(ffdata)

      feols(fml_base_w126, aadata,
        weights = ~fraction,
        nthreads = 0, notes = F, lean = T
      ) %>%
        tidy()
    }, .id = "id", .progress = T)
  }))

saveRDS(data %>% select(-fdata), "data/boots_f1_w126.rds")
