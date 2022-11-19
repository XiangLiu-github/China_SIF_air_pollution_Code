source("script/loadPackages.R")
source("script/loadFunctions.R")

mask <- rast("data/outputs/masks/mask_extraction.tif") %>% sum(na.rm = T)

a <- list.files("../data_archive/China_500m_irrigation/IrriMap_CN/", full.names = T) %>%
  pro_map(~ rast(.x) %>% resample(mask, method = "sum")) %>%
  rast()

b <- rast("../data_archive/China_500m_irrigation/IrriMap_CN/2000.tif")
b[is.na(b)] <- 1
b <- b %>% resample(mask, method = "sum")

a <- a / b
a[is.na(a)] <- 0

a <- a %>%
  as.data.frame(xy = T) %>%
  pivot_longer(-c(x, y),
    names_to = "year", names_transform = list(year = as.integer),
    values_to = "irg_fraction"
  ) %>%
  trim_xy()

check_join(a)

saveRDS(a, "data/outputs/irrigation/tidied.rds")
