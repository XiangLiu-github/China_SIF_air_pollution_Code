source("script/loadPackages.R")
source("script/loadFunctions.R")

mask <- rast("data/outputs/masks/mask_extraction.tif") %>% sum(na.rm = T)

data <- list.files("../data_archive/China_10km_ozone/Zhouetal2022/", full.names = T) %>%
  pro_map(function(afile) {
    r <- rast(afile) %>%
      `names<-`(time(.)) %>%
      flip()

    resample(r, mask, method = "near") %>%
      mask(mask) %>%
      as.data.frame(xy = T) %>%
      pivot_longer(-c(x, y),
        names_to = "date", names_transform = list(date = ymd),
        values_to = "O3"
      )
  })

data <- data %>%
  bind_rows() %>%
  trim_xy()

check_join(data)

require(parallel)
cl <- makeCluster(qn)

data_station <-
  inner_join(
    read_rds("../Vegetation_ozone/data/inputs/station_raw/ozone/China.rds"),
    read_rds("../pre_NS2021/data/inputs/station_raw/ozone/China.rds")
  ) %>%
  mutate(O3 = rollingo3 * 1e3 * 2)

# ug m-3 = ppb * 2

model_w126 <- bam(W126 ~ s(O3), data = data_station, family = tw(), chunk.size = 5000, cluster = cl)
model_aot40 <- bam(AOT40 ~ s(O3), data = data_station, family = tw(), chunk.size = 5000, cluster = cl)

data <-
  data %>%
  mutate(
    W126 = exp(predict(model_w126, ., cluster = cl)),
    AOT40 = exp(predict(model_aot40, ., cluster = cl)),
    year = year(date), month = month(date), .keep = "unused"
  )

saveRDS(data, "data/outputs/ozone/tidied.rds")
saveRDS(lst(model_aot40, model_w126), "data/outputs/ozone/gam.rds")
