source("script/loadPackages.R")
source("script/loadFunctions.R")

mask <- rast("data/outputs/masks/mask_extraction.tif") %>% sum(na.rm = T)

tidy_raster <- function(afile) {
  # the nc file to be read
  ncdata <- nc_open(afile)

  # read longitude and latitude from NetCDF file
  lon <- ncvar_get(ncdata, "longitude")
  lat <- ncvar_get(ncdata, "latitude")

  # read var data
  vars_1 <- ncvar_get(ncdata, "Cloud_Optical_Thickness_Total/Mean")
  vars_1[is.na(vars_1)] <- 0
  vars_2 <- ncvar_get(ncdata, "Cloud_Mask_Fraction/Mean")

  vars <- vars_1 * vars_2

  # close nc file
  nc_close(ncdata)

  ## (b) transform
  # convert var to raster in memory
  r <- rast(vars,
    crs = "+proj=longlat +datum=WGS84 +no_defs"
  )
  ext(r) <- c(min(lon) - 0.5, max(lon) + 0.5, min(lat) - 0.5, max(lat) + 0.5)

  r <- flip(r, direction = "vertical")

  names(r) <- "cloud"

  resample(r, mask, method = "near") %>%
    mask(mask) %>%
    as.data.frame(xy = T)
}

MODIS <-
  tibble(
    files = list.files("../data_archive/Global_100km_cloud/", ".nc", full.names = T),
    date = files %>% str_sub(48, 54) %>% parse_date_time("Yj")
  ) %>%
  mutate(
    data = pro_map(files, tidy_raster),
    year = year(date), month = month(date), .keep = "unused"
  )

MODIS <-
  MODIS %>%
  unnest() %>%
  trim_xy()

check_join(MODIS)

saveRDS(MODIS, "data/outputs/cloud/tidied.rds")
