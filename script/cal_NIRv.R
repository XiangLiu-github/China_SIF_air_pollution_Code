source("script/loadPackages.R")
source("script/loadFunctions.R")

# 1 km --------------------------------------------------------------------

mask <- rast("data/outputs/masks/mask_extraction_1km.tif") %>% sum(na.rm = T)

walk(2001:2020, function(ayear) {
  print(ayear)

  NIRv <-
    tibble(
      files = list.files("/Users/xiangliu/Downloads/", "tif", full.names = T),
      files_ = list.files("/Users/xiangliu/Downloads/", "tif")
    ) %>%
    separate(files_, c(NA, "year", "month", NA), "_|\\.", convert = T) %>%
    filter(year == ayear) %>%
    mutate(data = pro_map(files, function(afile) {
      adata <- rast(afile)
      project(adata, mask, method = "near") %>%
        mask(mask) %>%
        `names<-`("NIRv") %>%
        as.data.frame(xy = T)
    }))

  NIRv %>%
    select(-files) %>%
    unnest() %>%
    qsave(str_c("data/outputs/NIRv/", ayear, ".qs"), nthreads = 3)
})


# 0.05 --------------------------------------------------------------------

mask <- rast("data/outputs/masks/mask_extraction.tif") %>% sum(na.rm = T)

NIRv <-
  tibble(
    files = list.files("../data_archive/China_1km_NIRv/", full.names = T),
    files_ = list.files("../data_archive/China_1km_NIRv/"),
    date = str_remove(files_, "NIRv_") %>% str_remove(".tif") %>% ym()
  ) %>%
  select(-files_) %>%
  filter(year(date) %in% 2005:2019) %>%
  mutate(data = map(files, function(afile) {
    adata <- rast(afile)
    resample(adata, mask) %>%
      mask(mask) %>%
      `names<-`("NIRv") %>%
      as.data.frame(xy = T)
  }, .progress = T), .keep = "unused")

NIRv <- NIRv %>%
  mutate(year = year(date), month = month(date), .keep = "unused") %>%
  unnest() %>%
  trim_xy()

check_join(NIRv)

saveRDS(NIRv, "data/outputs/NIRv/tidied.rds")


# 0.05 new ----------------------------------------------------------------

mask <- rast("data/outputs/masks/mask_extraction.tif") %>% sum(na.rm = T)

rast("D:/Data/Global_0.05_Vegetation/MOD13C2.A2000032.061.2020042083108.hdf") %>%
  names()

MODIS <-
  tibble(
    files = list.files("D:/Data/Global_0.05_Vegetation/", full.names = T),
    files_ = list.files("D:/Data/Global_0.05_Vegetation/")
  ) %>%
  separate(files_, c(NA, "date", NA, NA, NA), "\\.") %>%
  mutate(date = str_remove(date, "A") %>% parse_date_time("Yj")) %>%
  mutate(
    year = year(date), month = month(date),
    data = map(files, function(afile) {
      a <- rast(afile)[[c("\"CMG 0.05 Deg Monthly NDVI\"", "\"CMG 0.05 Deg Monthly EVI\"", "\"CMG 0.05 Deg Monthly NIR reflectance\"")]] %>%
        `crs<-`(crs(mask)) %>%
        `names<-`(c("NDVI", "EVI", "NIR")) %>%
        crop(mask) %>%
        mask(mask)

      a <- a / 10000 / 10000

      a[["NIRv"]] <- a[["NIR"]] * a[["NDVI"]]
      a[["NIRvc"]] <- a[["NIR"]] * (a[["NDVI"]] - 0.08)


      a <- a[[c("NDVI", "EVI", "NIRv", "NIRvc")]]

      a %>% as.data.frame(xy = T)
    }, .progress = T), .keep = "unused"
  ) %>%
  unnest() %>%
  trim_xy()

saveRDS(MODIS, "data/outputs/NIRv/tidied.rds")
