source("script/loadPackages.R")
source("script/loadFunctions.R")

mask <- rast("data/outputs/masks/mask_extraction.tif") %>% sum(na.rm = T)

############################################
# mW m−2 nm−1 sr−1 monthly
GOSIF <-
  expand_grid(
    year = 2000:2021,
    month = str_pad(1:12, width = 2, side = "left", pad = "0")
  ) %>%
  filter(!(year == 2000 & month %in% c("01", "02"))) %>%
  mutate(data = pro_map2(year, month, function(ayear, amonth) {
    afile <- str_c("/vsigzip/../data_archive/Global_5km_SIF/GOSIF_monthly/", ayear, "_", amonth, ".tif.gz")

    arst <- rast(afile) %>%
      `ext<-`(c(-180, 180, -90, 90)) %>%
      resample(mask, method = "near") %>%
      mask(mask) %>%
      `names<-`("GOSIF")

    arst[arst >= 32766] <- NA

    as.data.frame(arst, xy = T)
  }))

GOSIF <- GOSIF %>%
  mutate(month = as.integer(month)) %>%
  unnest() %>%
  trim_xy()

check_join(GOSIF)

saveRDS(GOSIF, "data/outputs/SIF/GOSIF.rds")

############################################
# mW m−2 nm−1 sr−1 8 day

GOSIF <- tibble(
  files_ = list.files("../data_archive/Global_5km_SIF/GOSIF_8day/"),
  files = list.files("../data_archive/Global_5km_SIF/GOSIF_8day/", full.names = T)
) %>%
  mutate(
    date = str_remove(files_, "GOSIF_") %>%
      str_remove("\\.tif\\.gz") %>%
      parse_date_time("Yj"),
    files = str_c("/vsigzip/", files), .keep = "unused"
  ) %>%
  mutate(data = pro_map(files, function(afile) {
    arst <- rast(afile) %>%
      `ext<-`(c(-180, 180, -90, 90)) %>%
      crop(mask) %>%
      mask(mask) %>%
      `names<-`("GOSIF")

    arst[arst >= 32766] <- NA

    as.data.frame(arst, xy = T)
  }))

GOSIF <- GOSIF %>%
  select(-files) %>%
  unnest()

saveRDS(GOSIF, "data/outputs/SIF/GOSIF_8day.rds")

############################################
# mW m−2 nm−1 sr−1

RTSIF <-
  tibble(
    files_ = list.files("../data_archive/Global_5km_SIF/RTSIF_8day/"),
    files = list.files("../data_archive/Global_5km_SIF/RTSIF_8day/", full.names = T)
  ) %>%
  separate(files_, c(NA, "Date", NA), "_|\\.") %>%
  mutate(
    Date = ymd(Date),
    YM = str_c(year(Date), "_", month(Date)),
    data = map(files, rast)
  )

RTSIF_data <- rast(RTSIF$data)

RTSIF <- tapp(RTSIF_data, RTSIF$YM, mean)

RTSIF <-
  RTSIF %>%
  resample(mask, method = "near") %>%
  mask(mask)

RTSIF <-
  RTSIF %>%
  as.data.frame(xy = T) %>%
  pivot_longer(-c(x, y),
    names_prefix = "X", names_to = c("year", "month"), names_sep = "_",
    names_transform = list(year = as.integer, month = as.integer),
    values_to = "RTSIF"
  ) %>%
  trim_xy()

check_join(RTSIF)

saveRDS(RTSIF, "data/outputs/SIF/RTSIF.rds")

############################################
# mW m−2 nm−1 sr−1

CSIF <-
  tibble(
    files_ = list.files("../data_archive/Global_5km_SIF/CSIF_monthly//", recursive = T),
    files = list.files("../data_archive/Global_5km_SIF/CSIF_monthly//", recursive = T, full.names = T)
  ) %>%
  separate(files_, c(NA, NA, NA, NA, "date", NA, NA), "\\.") %>%
  mutate(date = parse_date_time(date, "Yj")) %>%
  filter(date >= ymd("2000-03-01")) %>%
  mutate(
    year = year(date), month = month(date),
    data = pro_map(files, function(afile) {
      temp <- raster(afile, varname = "clear_daily_SIF")[]

      rast(resolution = c(0.05, 0.05)) %>%
        terra::`values<-`(temp) %>%
        resample(mask, method = "near") %>%
        mask(mask)
    }), .keep = "unused"
  )


CSIF_data <- tapp(rast(CSIF$data), str_c(CSIF$year, "_", CSIF$month), mean)

CSIF_data <- CSIF_data %>%
  as.data.frame(xy = T) %>%
  pivot_longer(-c(x, y),
    names_prefix = "X", names_to = c("year", "month"), names_sep = "_",
    names_transform = list(year = as.integer, month = as.integer),
    values_to = "CSIF"
  ) %>%
  trim_xy()

check_join(CSIF_data)

saveRDS(CSIF_data, "data/outputs/SIF/CSIF.rds")
