source("script/loadPackages.R")
source("script/loadFunctions.R")

exts <- rast("../data_archive/China_1km_majorcrop_area/CHN_Maize_2019.tif") %>% ext()

# make a unified 0.05 mask for extraction  --------------------------------

all <-
  list.files("../data_archive/China_1km_majorcrop_area", full.names = T) %>%
  map(rast) %>%
  pro_map(~ extend(.x, exts)) %>%
  rast() %>%
  tapp(rep(c("Maize", "Rice", "Wheat"), each = 20), "sum", na.rm = T)

all[all > 0] <- 1
all[all == 0] <- NA

all <- all %>%
  project("epsg:4326", method = "near")

degree <- rast(resolution = 0.05) %>%
  crop(all, snap = "out")

all <- all %>% resample(degree, "sum")

writeRaster(all, "data/outputs/masks/mask_extraction.tif", overwrite = TRUE)


# make annual 0.05 mask for joining ---------------------------------------

c("Maize", "Wheat", "Rice") %>%
  walk(function(acrop) {
    print(acrop)

    all <-
      list.files("../data_archive/China_1km_majorcrop_area", acrop, full.names = T) %>%
      map(rast) %>%
      map(~ extend(.x, exts)) %>%
      rast() %>%
      `names<-`(2000:2019) %>%
      project("epsg:4326", method = "near")

    degree <- rast(resolution = 0.05) %>%
      crop(all, snap = "out")

    all <- all %>% resample(degree, "sum")

    all[all == 0] <- NA

    writeRaster(all, str_c("data/outputs/masks/mask_", acrop, ".tif"), overwrite = T)
  })

# make a unified 1km mask for extraction  --------------------------------

temp_1km <- rast(resolution = 0.01) %>%
  crop(rast("data/outputs/masks/mask_extraction.tif"), snap = "out")

all <-
  list.files("../data_archive/China_1km_majorcrop_area", full.names = T) %>%
  map(rast) %>%
  pro_map(~ extend(.x, exts)) %>%
  rast() %>%
  tapp(rep(c("Maize", "Rice", "Wheat"), each = 20), "sum", na.rm = T)

all[all > 0] <- 1
all[all == 0] <- NA

all <- project(all, temp_1km, method = "near")

writeRaster(all, "data/outputs/masks/mask_extraction_1km.tif", overwrite = TRUE)
