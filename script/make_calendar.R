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
