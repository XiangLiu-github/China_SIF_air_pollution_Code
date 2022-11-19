source("script/loadPackages.R")
source("script/loadFunctions.R")

mask <- rast("data/outputs/masks/mask_extraction.tif") %>%
  sum(na.rm = T)

# temperature bins --------------------------------------------------------

# qian wan zhi neng zai wan shang yun xing ~ 8 h

temp_bins <-
  tibble(
    files_ = list.files("D:/Data/China_1km_climate/temp_bins/", ".tif"),
    files = list.files("D:/Data/China_1km_climate/temp_bins/", ".tif", full.names = T)
  ) %>%
  separate(files_, c("year", "month", NA), "_|\\.", convert = T) %>%
  arrange(year, month) %>%
  mutate(data = pro_map(files, function(afile) {
    rast(afile) %>%
      resample(mask, method = "average") %>%
      mask(mask) %>%
      as.data.frame(xy = T) %>%
      trim_xy()
  }), .keep = "unused")

# temp_bins = temp_bins %>% unnest()

check_join(temp_bins$data[[5]])

qsave(temp_bins, "data/outputs/temp/temp_bins.qs", nthreads = 5)


# tmax mean ---------------------------------------------------------------
tmax <-
  stack("../data_archive/China_1km_climate/China_1km_maxtmp_monmean.nc") %>%
  resample(mask %>% raster()) %>%
  mask(mask %>% raster()) %>%
  as.data.frame(xy = T) %>%
  trim_xy()

tmax <- tmax %>%
  pivot_longer(-c(x, y),
    names_to = "date", names_prefix = "X", names_transform = list(date = ymd),
    values_to = "maxtmp", values_drop_na = T
  ) %>%
  mutate(year = year(date), month = month(date), .keep = "unused")

check_join(tmax)

saveRDS(tmax, "data/outputs/temp/tmax.rds")

