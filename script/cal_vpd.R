source("script/loadPackages.R")
source("script/loadFunctions.R")

mask <- rast("data/outputs/masks/mask_extraction.tif") %>% sum(na.rm = T)

VPD <- tibble(
  year = 2005:2019,
  data = map(year, function(ayear) {
    P <- (rast(str_c("../data_archive/ERA5_0.25_climate/e5.moda.an.sfc.128_151_msl.ll025sc.", ayear, "010100_", ayear, "120100.nc")) / 100) %>% # Pa -> hPa
      `crs<-`("+proj=longlat +datum=WGS84") %>%
      project(mask) %>%
      mask(mask)
    Tt <- (rast(str_c("../data_archive/ERA5_0.25_climate/e5.moda.an.sfc.128_167_2t.ll025sc.", ayear, "010100_", ayear, "120100.nc")) - 273.15) %>% # K -> Cdegree
      `crs<-`("+proj=longlat +datum=WGS84") %>%
      project(mask) %>%
      mask(mask)
    Td <- (rast(str_c("../data_archive/ERA5_0.25_climate/e5.moda.an.sfc.128_168_2d.ll025sc.", ayear, "010100_", ayear, "120100.nc")) - 273.15) %>% # K -> Cdegree
      `crs<-`("+proj=longlat +datum=WGS84") %>%
      project(mask) %>%
      mask(mask)

    # ref: https://www.science.org/doi/10.1126/sciadv.aax1396; https://doi.org/10.1029/2022EF003019

    fw <- 1 + 7e-4 + 3.46 * 1e-6 * P

    SVP <- 6.112 * fw * exp(17.67 * Tt / (Tt + 243.5))
    AVP <- 6.112 * fw * exp(17.67 * Td / (Td + 243.5))
    VPD <- SVP - AVP
    names(VPD) <- 1:12

    VPD %>% as.data.frame(xy = T)
  }, .progress = T)
)

VPD <- VPD %>%
  unnest() %>%
  pivot_longer(-c(x, y, year),
    names_to = "month", names_transform = list(month = as.integer),
    values_to = "vpd"
  ) %>%
  trim_xy()

check_join(VPD)

saveRDS(VPD, "data/outputs/vpd/tidied.rds")


a <- qread("../tidy_crop_statistics/data/tidied.qs") %>%
  st_drop_geometry() %>%
  filter(crop == "玉米", year %in% 2005:2019)
