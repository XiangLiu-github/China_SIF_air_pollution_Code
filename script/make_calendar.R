source("script/loadPackages.R")
source("script/loadFunctions.R")

# make a unified 0.05 calendar for joining --------------------------------
# monthly
exts <- rast("../data_archive/China_1km_majorcrop_phenology/CHN_Maize_HE_2019.tif") %>% ext()

datap <-
  tibble(
    period_ = c(
      "Maize_HE", "Maize_MA", "Maize_V3",
      "Rice\\(LR\\)_HE", "Rice\\(LR\\)_MA", "Rice\\(LR\\)_TR",
      "Rice\\(SR&ER\\)_HE", "Rice\\(SR&ER\\)_MA", "Rice\\(SR&ER\\)_TR",
      "Wheat_GR&EM", "Wheat_HE", "Wheat_MA"
    ),
    period = c(
      "Maize_HE", "Maize_MA", "Maize_V3",
      "Rice(LR)_HE", "Rice(LR)_MA", "Rice(LR)_TR",
      "Rice(SR&ER)_HE", "Rice(SR&ER)_MA", "Rice(SR&ER)_TR",
      "Wheat_GR&EM", "Wheat_HE", "Wheat_MA"
    )
  ) %>%
  mutate(
    period = str_replace_all(period, "V3|TR", "GR&EM"),
    data = map2(period_, period, function(aperiod, ap) {
      print(aperiod)

      all <- list.files("../data_archive/China_1km_majorcrop_phenology/", aperiod, full.names = T) %>%
        map(rast) %>%
        map(~ extend(.x, exts)) %>%
        rast() %>%
        `names<-`(2000:2019) %>%
        project("epsg:4326", method = "near")

      degree <- rast(resolution = 0.05) %>%
        crop(all, snap = "out")

      all <- all %>% resample(degree, "average")

      all <- all %>%
        as.data.frame(xy = T, na.rm = F) %>%
        pivot_longer(-c(x, y),
          values_drop_na = T,
          names_to = "year", names_transform = list(year = as.integer)
        ) %>%
        mutate(
          value = round(value) %>%
            as.character() %>%
            str_c(year, "-", .) %>%
            parse_date_time("Y-j"),
          month = month(value), day = day(value)
        )

      if (str_detect(ap, "GR&EM")) {
        all <- all %>%
          mutate(month_fnl = fifelse(day >= 15, month + 1, month))
      } else if (str_detect(ap, "MA")) {
        all <- all %>%
          mutate(month_fnl = fifelse(day >= 15, month, month - 1))
      } else if (str_detect(ap, "HE")) {
        all <- all %>%
          mutate(month_fnl = month)
      }

      all %>% select(x, y, year, month_fnl)
    })
  ) %>%
  select(-period_)

datap <-
  datap %>%
  separate(period, c("crop", "stage"), "_") %>%
  unnest() %>%
  pivot_wider(names_from = "stage", values_from = month_fnl) %>%
  mutate(month = map2(`GR&EM`, MA, seq)) %>%
  unnest()

datap <-
  datap %>%
  filter(MA >= `GR&EM`)

saveRDS(datap, "data/outputs/calendar/tidied.rds")


# 8day --------------------------------------------------------------------

exts <- rast("../data_archive/China_1km_majorcrop_phenology/CHN_Maize_HE_2019.tif") %>% ext()

days_tar <- list.files("../data_archive/Global_5km_SIF/GOSIF_8day/", "2020") %>%
  str_remove("GOSIF_2020") %>%
  str_remove(".tif.gz") %>%
  as.integer()

cuts <- c(-Inf, days_tar[-1] - diff(days_tar) / 2, Inf)

datap <-
  tibble(
    period_ = c(
      "Maize_HE", "Maize_MA", "Maize_V3",
      "Rice\\(LR\\)_HE", "Rice\\(LR\\)_MA", "Rice\\(LR\\)_TR",
      "Rice\\(SR&ER\\)_HE", "Rice\\(SR&ER\\)_MA", "Rice\\(SR&ER\\)_TR",
      "Wheat_GR&EM", "Wheat_HE", "Wheat_MA"
    ),
    period = c(
      "Maize_HE", "Maize_MA", "Maize_V3",
      "Rice(LR)_HE", "Rice(LR)_MA", "Rice(LR)_TR",
      "Rice(SR&ER)_HE", "Rice(SR&ER)_MA", "Rice(SR&ER)_TR",
      "Wheat_GR&EM", "Wheat_HE", "Wheat_MA"
    )
  ) %>%
  mutate(data = map(period_, function(aperiod) {
    print(aperiod)

    temp <- raster("../data_archive/China_1km_climate/China_1km_maxtmp_monmean.nc")
    temp <- rast(
      nrows = dim(temp)[1], ncols = dim(temp)[2],
      xmin = extent(temp)@xmin, xmax = extent(temp)@xmax,
      ymin = extent(temp)@ymin, ymax = extent(temp)@ymax,
      crs = crs(temp)
    )

    all <- list.files("../data_archive/China_1km_majorcrop_phenology/", aperiod, full.names = T) %>%
      map(rast) %>%
      map(~ extend(.x, exts)) %>%
      rast() %>%
      `names<-`(2000:2019) %>%
      project(temp,
        method = "near"
      )

    degree <- rast(resolution = 0.05) %>%
      crop(all, snap = "out")

    all <- all %>% resample(degree, "average")

    all %>%
      as.data.frame(xy = T, na.rm = F) %>%
      pivot_longer(-c(x, y),
        values_drop_na = T,
        names_to = "year", names_transform = list(year = as.integer)
      ) %>%
      mutate(
        doy = days_tar[findInterval(value, cuts)], # find nearest value https://stackoverflow.com/questions/43472234/fastest-way-to-find-nearest-value-in-vector
      ) %>%
      select(-value)
  }), .keep = "unused")

test <- datap %>%
  mutate(period = str_replace_all(period, "V3|TR", "GR&EM")) %>%
  separate(period, c("crop", "stage"), "_") %>%
  unnest() %>%
  mutate(date = parse_date_time(str_c("0001_", doy), "Y_j") + years(year - 1)) %>%
  select(-c(doy))

test <- test %>%
  pivot_wider(names_from = "stage", values_from = date) %>%
  filter(MA > `GR&EM`) %>%
  mutate(date = map2(`GR&EM`, MA, ~ seq(.x, .y, by = "8 days"))) %>%
  unnest()

saveRDS(test, "data/outputs/calendar/tidied_8day.rds")







# make a unified 1 km calendar for joining --------------------------------

mask <- rast("data/outputs/masks/mask_extraction_1km.tif") %>%
  sum(na.rm = T) %>%
  wrap()

pro_walk(2000:2019, function(ayear) {
  plan(multisession, workers = 4)

  datap <-
    tibble(
      period_ = c(
        "Maize_HE", "Maize_MA", "Maize_V3",
        "Rice\\(LR\\)_HE", "Rice\\(LR\\)_MA", "Rice\\(LR\\)_TR",
        "Rice\\(SR&ER\\)_HE", "Rice\\(SR&ER\\)_MA", "Rice\\(SR&ER\\)_TR",
        "Wheat_GR&EM", "Wheat_HE", "Wheat_MA"
      ),
      period = c(
        "Maize_HE", "Maize_MA", "Maize_V3",
        "Rice(LR)_HE", "Rice(LR)_MA", "Rice(LR)_TR",
        "Rice(SR&ER)_HE", "Rice(SR&ER)_MA", "Rice(SR&ER)_TR",
        "Wheat_GR&EM", "Wheat_HE", "Wheat_MA"
      )
    ) %>%
    mutate(data = future_map(period_, function(aperiod) {
      all <- list.files("../data_archive/China_1km_majorcrop_phenology/", str_c(aperiod, "_", ayear), full.names = T) %>%
        rast() %>%
        `names<-`(ayear) %>%
        project(rast(mask), method = "near")

      all %>%
        as.data.frame(xy = T) %>%
        pivot_longer(-c(x, y),
          names_to = "year", names_transform = list(year = as.integer)
        ) %>%
        mutate(
          value = round(value) %>%
            as.character() %>%
            str_c(year, "-", .) %>%
            parse_date_time("Y-j"),
          month = month(value), day = day(value),
          month_fnl = fifelse(day >= 15, month + 1, month)
        ) %>%
        select(x, y, year, month_fnl)
    }, .progress = T), .keep = "unused")

  plan(sequential)

  datap <-
    datap %>%
    mutate(period = str_replace_all(period, "V3|TR", "GR&EM")) %>%
    separate(period, c("crop", "stage"), "_") %>%
    unnest() %>%
    pivot_wider(names_from = "stage", values_from = month_fnl) %>%
    mutate(month = map2(`GR&EM`, MA, seq)) %>%
    unnest()

  datap <-
    datap %>%
    filter(month <= 12)

  qsave(datap, str_c("data/outputs/calendar/tidied_1km_", ayear, ".qs"), nthreads = 5)
})
