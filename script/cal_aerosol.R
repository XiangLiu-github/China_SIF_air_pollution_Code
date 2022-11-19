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

