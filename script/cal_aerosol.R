source("script/loadPackages.R")
source("script/loadFunctions.R")

mask <- rast("data/outputs/masks/mask_extraction.tif") %>%
  sum(na.rm = T) %>%
  wrap()

# AOD ---------------------------------------------------------------------

plan(multisession, workers = 3)

AOD <-
  tibble(
    files_ = list.files("D:/Data/China_1km_air_pollution/LGHAP.AOD.monthly/"),
    files = list.files("D:/Data/China_1km_air_pollution/LGHAP.AOD.monthly/", full.names = T)
  ) %>%
  separate(files_, c(NA, NA, NA, "date", NA), "\\.") %>%
  mutate(date = str_remove(date, "M") %>% ym()) %>%
  mutate(year = year(date), month = month(date), .keep = "unused") %>%
  mutate(data = future_map(files, function(afile) {
    ncdata <- nc_open(afile)

    lon <- ncvar_get(ncdata, "lon")
    lat <- ncvar_get(ncdata, "lat")

    vars <- ncvar_get(ncdata, "AOD")

    nc_close(ncdata)

    r <- rast(vars,
      crs = "+proj=longlat +datum=WGS84 +no_defs"
    )
    ext(r) <- c(min(lon), max(lon), min(lat), max(lat))

    r <- flip(r, direction = "vertical")

    names(r) <- "AOD"

    resample(r, rast(mask), method = "average") %>%
      mask(rast(mask)) %>%
      as.data.frame(xy = T)
  },
  .progress = T,
  .options = furrr_options(packages = c("ncdf4", "terra", "tidyverse"))
  ), .keep = "unused")

plan(sequential)

AOD <- AOD %>%
  unnest() %>%
  trim_xy()

check_join(AOD)

saveRDS(AOD, "data/outputs/aerosol/AOD.rds")

# PM25 --------------------------------------------------------------------

plan(multisession, workers = 3)

PM25 <-
  tibble(
    files_ = list.files("D:/Data/China_1km_air_pollution/LGHAP.PM25.monthly/"),
    files = list.files("D:/Data/China_1km_air_pollution/LGHAP.PM25.monthly/", full.names = T)
  ) %>%
  separate(files_, c(NA, NA, NA, "date", NA), "\\.") %>%
  mutate(date = str_remove(date, "M") %>% ym()) %>%
  mutate(year = year(date), month = month(date), .keep = "unused") %>%
  mutate(data = future_map(files, function(afile) {
    ncdata <- nc_open(afile)

    lon <- ncvar_get(ncdata, "lon")
    lat <- ncvar_get(ncdata, "lat")

    vars <- ncvar_get(ncdata, "PM25")

    nc_close(ncdata)

    r <- rast(vars,
      crs = "+proj=longlat +datum=WGS84 +no_defs"
    )
    ext(r) <- c(min(lon), max(lon), min(lat), max(lat))

    r <- flip(r, direction = "vertical")

    names(r) <- "PM25"

    resample(r, rast(mask), method = "average") %>%
      mask(rast(mask)) %>%
      as.data.frame(xy = T)
  },
  .progress = T,
  .options = furrr_options(packages = c("ncdf4", "terra", "tidyverse"))
  ), .keep = "unused")

plan(sequential)

PM25 <- PM25 %>%
  unnest() %>%
  trim_xy()

check_join(PM25)

saveRDS(PM25, "data/outputs/aerosol/PM25.rds")


# PM10 --------------------------------------------------------------------

plan(multisession, workers = 3)

PM10 <-
  tibble(
    files_ = list.files("D:/Data/China_1km_air_pollution/LGHAP.PM10.monthly/"),
    files = list.files("D:/Data/China_1km_air_pollution/LGHAP.PM10.monthly/", full.names = T)
  ) %>%
  separate(files_, c(NA, NA, NA, "date", NA), "\\.") %>%
  mutate(date = str_remove(date, "M") %>% ym()) %>%
  mutate(year = year(date), month = month(date), .keep = "unused") %>%
  mutate(data = future_map(files, function(afile) {
    ncdata <- nc_open(afile)

    lon <- ncvar_get(ncdata, "lon")
    lat <- ncvar_get(ncdata, "lat")

    vars <- ncvar_get(ncdata, "PM10")

    nc_close(ncdata)

    r <- rast(vars,
      crs = "+proj=longlat +datum=WGS84 +no_defs"
    )
    ext(r) <- c(min(lon), max(lon), min(lat), max(lat))

    r <- flip(r, direction = "vertical")

    names(r) <- "PM10"

    resample(r, rast(mask), method = "average") %>%
      mask(rast(mask)) %>%
      as.data.frame(xy = T)
  },
  .progress = T,
  .options = furrr_options(packages = c("ncdf4", "terra", "tidyverse"))
  ), .keep = "unused")

plan(sequential)

PM10 <- PM10 %>%
  unnest() %>%
  trim_xy()

check_join(PM10)

saveRDS(PM10, "data/outputs/aerosol/PM10.rds")

# AOD 1km ---------------------------------------------------------------------

mask <- rast("data/outputs/masks/mask_extraction_1km.tif") %>%
  sum(na.rm = T) %>%
  wrap()

pro_walk(2000:2020, function(ayear) {
  plan(multisession, workers = 4)

  AOD <-
    tibble(
      files_ = list.files("D:/Data/China_1km_air_pollution/LGHAP.AOD.monthly/"),
      files = list.files("D:/Data/China_1km_air_pollution/LGHAP.AOD.monthly/", full.names = T)
    ) %>%
    separate(files_, c(NA, NA, NA, "date", NA), "\\.") %>%
    mutate(date = str_remove(date, "M") %>% ym()) %>%
    mutate(year = year(date), month = month(date), .keep = "unused") %>%
    filter(year == ayear) %>%
    mutate(data = future_map(files, function(afile) {
      ncdata <- nc_open(afile)

      lon <- ncvar_get(ncdata, "lon")
      lat <- ncvar_get(ncdata, "lat")

      vars <- ncvar_get(ncdata, "AOD")

      nc_close(ncdata)

      r <- rast(vars,
        crs = "+proj=longlat +datum=WGS84 +no_defs"
      )
      ext(r) <- c(min(lon), max(lon), min(lat), max(lat))

      r <- flip(r, direction = "vertical")

      names(r) <- "AOD"

      project(r, rast(mask), method = "near") %>%
        mask(rast(mask)) %>%
        as.data.frame(xy = T)
    },
    .progress = T,
    .options = furrr_options(packages = c("ncdf4", "terra", "tidyverse"))
    ), .keep = "unused")

  plan(sequential)

  AOD <- AOD %>% unnest()

  qsave(AOD, str_c("data/outputs/aerosol/AOD_1km_", ayear, ".qs"), nthreads = 5)
})


# MAIAC -------------------------------------------------------------------


plan(multisession, workers = 3)

MAIAC <- tibble(
  files_ = list.files("D:/Data/China_1km_air_pollution/MAIAC_AOD_monthly/"),
  files = list.files("D:/Data/China_1km_air_pollution/MAIAC_AOD_monthly/", full.names = T)
) %>%
  separate(files_, c(NA, NA, "date", NA), "_|\\.") %>%
  mutate(date = date %>% ym()) %>%
  mutate(year = year(date), month = month(date), .keep = "unused") %>%
  mutate(data = future_map(files, function(afile) {
    r <- rast(afile)

    r[r == 0] <- NA

    names(r) <- "MAIAC"

    resample(r, rast(mask), method = "average") %>%
      mask(rast(mask)) %>%
      as.data.frame(xy = T)
  },
  .progress = T,
  .options = furrr_options(packages = c("terra", "tidyverse"))
  ), .keep = "unused")

plan(sequential)

MAIAC <- MAIAC %>%
  unnest() %>%
  trim_xy()

check_join(MAIAC)

saveRDS(MAIAC, "data/outputs/aerosol/MAIAC.rds")
