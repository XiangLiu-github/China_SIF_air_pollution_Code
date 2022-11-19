source("script/loadPackages.R")
source("script/loadFunctions.R")

mask <- rast("data/outputs/masks/mask_extraction.tif") %>% sum(na.rm = T)

# mean to daily -----------------------------------------------------------

walk(2000:2020, function(ayear) {
  print(ayear)

  hourly_data <-
    tibble(
      files_ = list.files("D:/Data/China_25km_soilmoisture/", ".nc4"),
      files = list.files("D:/Data/China_25km_soilmoisture/", ".nc4", full.names = T)
    ) %>%
    separate(files_, c(NA, "date", "hour", NA, NA, NA, NA), "\\.") %>%
    mutate(date = str_remove(date, "A") %>% ymd()) %>%
    filter(year(date) == ayear) %>%
    mutate(
      data_s = map(files, ~ rast(.x)[["SoilMoi0_10cm_inst"]]),
      data_a = map(files, ~ rast(.x)[["RootMoist_inst"]])
    )

  surface <- rast(hourly_data$data_s) %>%
    tapp(hourly_data$date, "mean")

  root <- rast(hourly_data$data_a) %>%
    tapp(hourly_data$date, "mean")

  writeRaster(surface, str_c("data/inputs/soil_moisture/surface_", ayear, ".tif"), overwrite = T)
  writeRaster(root, str_c("data/inputs/soil_moisture/root_", ayear, ".tif"), overwrite = T)
})


# daily to bins -----------------------------------------------------------

sm_bins <-
  tibble(
    files_ = list.files("data/inputs/soil_moisture/"),
    files = list.files("data/inputs/soil_moisture/", full.names = T)
  ) %>%
  filter(!str_detect(files_, "aux")) %>%
  separate(files_, c("type", NA), "_") %>%
  mutate(data = pro_map2(files, type, function(afile, atype) {
    if (atype == "root") {
      bins <- c(seq(0, 800, by = 100), 10000)
    } else {
      bins <- c(seq(0, 40, by = 5), 1000)
    }

    rast(afile) %>%
      classify(bins, include.lowest = TRUE, brackets = TRUE) %>%
      as.data.frame(xy = T) %>%
      pivot_longer(-c(x, y),
        names_prefix = "X", names_to = "date",
        names_transform = list(date = ymd)
      ) %>%
      lazy_dt() %>%
      mutate(year = year(date), month = month(date), .keep = "unused") %>%
      count(x, y, value, year, month) %>%
      pivot_wider(names_from = value, values_from = n, values_fill = 0) %>%
      as_tibble()
  }), .keep = "unused")

root <-
  pro_map(sm_bins$data[1:21], function(adf) {
    adf %>%
      pivot_wider(
        names_from = c(year, month),
        values_from = c(
          `[0–100]`, `(100–200]`, `(200–300]`,
          `(300–400]`, `(400–500]`, `(500–600]`,
          `(600–700]`, `(700–800]`, `(800–10000]`
        )
      ) %>%
      rast(type = "xyz", crs = crs(mask)) %>%
      resample(mask, method = "near") %>%
      mask(mask) %>%
      as.data.frame(xy = T) %>%
      pivot_longer(-c(x, y),
        names_to = c("bins", "year", "month"), names_sep = "_",
        names_transform = list(year = as.integer, month = as.integer)
      ) %>%
      pivot_wider(names_from = bins) %>%
      trim_xy()
  })

root <- bind_rows(root) %>%
  `names<-`(c("x", "y", "year", "month", str_c("root_", 1:9)))

check_join(root)

saveRDS(root, "data/outputs/soil_moisture/root.rds")

surface <-
  pro_map(sm_bins$data[22:42], function(adf) {
    adf %>%
      pivot_wider(
        names_from = c(year, month),
        values_from = c(
          `[0–5]`, `(5–10]`, `(10–15]`,
          `(15–20]`, `(20–25]`, `(25–30]`,
          `(30–35]`, `(35–40]`, `(40–1000]`
        )
      ) %>%
      rast(type = "xyz", crs = crs(mask)) %>%
      resample(mask, method = "near") %>%
      mask(mask) %>%
      as.data.frame(xy = T) %>%
      pivot_longer(-c(x, y),
        names_to = c("bins", "year", "month"), names_sep = "_",
        names_transform = list(year = as.integer, month = as.integer)
      ) %>%
      pivot_wider(names_from = bins) %>%
      trim_xy()
  })

surface <- bind_rows(surface) %>%
  `names<-`(c("x", "y", "year", "month", str_c("surface_", 1:9)))

check_join(surface)

saveRDS(surface, "data/outputs/soil_moisture/surface.rds")


# month mean --------------------------------------------------------------


month_mean <-
  tibble(
    files_ = list.files("data/inputs/soil_moisture/"),
    files = list.files("data/inputs/soil_moisture/", full.names = T)
  ) %>%
  filter(!str_detect(files_, "aux")) %>%
  separate(files_, c("type", "year", NA), "_|\\.", convert = T) %>%
  mutate(data = pro_map(files, function(afile) {
    temp <- rast(afile)

    dates <- ymd(str_remove(names(temp), "X"))

    tapp(temp, str_c(year(dates), "_", month(dates)), mean) %>%
      resample(mask, method = "near") %>%
      mask(mask) %>%
      as.data.frame(xy = T) %>%
      pivot_longer(-c(x, y),
        names_to = c("year", "month"), names_prefix = "X", names_sep = "_", names_transform = list(year = as.integer, month = as.integer)
      )
  }))

month_mean <-
  month_mean %>%
  select(type, data) %>%
  unnest() %>%
  trim_xy() %>%
  pivot_wider(names_from = type)

check_join(month_mean)

saveRDS(month_mean, "data/outputs/soil_moisture/mean.rds")
