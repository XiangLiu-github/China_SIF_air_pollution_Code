source("script/loadPackages.R")
source("script/loadFunctions.R")

mask <- rast("data/outputs/masks/mask_extraction.tif") %>% sum(na.rm = T)

data <- list.files("../data_archive/China_10km_ozone/Zhouetal2022/", full.names = T) %>%
  map(function(afile) {
    r <- rast(afile) %>%
      `names<-`(time(.))

    resample(r, mask, method = "near") %>%
      mask(mask) %>%
      as.data.frame(xy = T) %>%
      pivot_longer(-c(x, y),
        names_to = "date", names_transform = list(date = ymd),
        values_to = "O3"
      )
  }, .progress = T)

data <- data %>%
  bind_rows() %>%
  filter(year(date) %in% 2005:2019) %>%
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

# 1km ---------------------------------------------------------------------

mask <- rast("data/outputs/masks/mask_extraction_1km.tif") %>%
  sum(na.rm = T)

require(parallel)
nc <- 9 ## cluster size, set for example portability
if (detectCores() > 1) { ## no point otherwise
  cl <- makeCluster(nc)
  ## could also use makeForkCluster, but read warnings first!
} else {
  cl <- NULL
}

pred_models <- read_rds("data/outputs/ozone/gam.rds")

list.files("D:/Data/China_10km_ozone/Zhouetal2022/", full.names = T) %>%
  pro_walk(function(afile) {
    r <- rast(afile) %>%
      `names<-`(time(.)) %>%
      flip()

    O3 <-
      resample(r, mask, method = "near") %>%
      mask(mask) %>%
      as.data.frame(xy = T) %>%
      pivot_longer(-c(x, y),
        names_to = "date", names_transform = list(date = ymd),
        values_to = "O3"
      ) %>%
      mutate(
        W126 = exp(predict(pred_models$model_w126, ., cluster = cl)),
        AOT40 = exp(predict(pred_models$model_aot40, ., cluster = cl)),
        year = year(date), month = month(date), .keep = "unused"
      )

    qsave(O3, str_c("data/outputs/ozone/tidied_", unique(O3$year), ".qs"), nthreads = 5)
  })
