source("script/loadPackages.R")
source("script/loadFunctions.R")

mask <- rast("data/outputs/masks/mask_extraction.tif") %>% sum(na.rm = T)

prep <-
  rast("../data_archive/China_1km_climate/China_1km_prep_monsum.nc") %>%
  `names<-`(time(.)) %>%
  resample(mask, method = "average") %>%
  mask(mask) %>%
  as.data.frame(xy = T) %>%
  trim_xy()

prep <- prep %>%
  pivot_longer(-c(x, y),
    names_to = "date", names_transform = list(date = ymd),
    values_to = "prep", values_drop_na = T
  ) %>%
  mutate(year = year(date), month = month(date), .keep = "unused")

check_join(prep)

saveRDS(prep, "data/outputs/prep/tidied.rds")
