source("script/loadPackages.R")
source("script/loadFunctions.R")

require(parallel)
cl <- makeCluster(qn)

cldr <-
  read_rds("data/outputs/calendar/tidied.rds") %>%
  filter((MA - `GR&EM`) >= 2) %>%
  filter((MA - month) <= 2) %>%
  trim_xy()

fraction <-
  tibble(
    crop = c("Maize", "Rice(LR)", "Rice(SR&ER)", "Wheat"),
    data = map(c("Maize", "Rice", "Rice", "Wheat"), function(acrop) {
      afile <- str_c("data/outputs/masks/mask_", acrop, ".tif")
      rast(afile) %>%
        as.data.frame(xy = T, na.rm = F) %>%
        pivot_longer(-c(x, y),
          names_to = "year", names_transform = list(year = as.integer),
          values_drop_na = T, values_to = "fraction"
        )
    })
  ) %>%
  unnest() %>%
  trim_xy()

ozone <- read_rds("data/outputs/ozone/tidied.rds")

ozone <- cldr %>%
  lazy_dt() %>%
  inner_join(ozone) %>%
  group_by(crop, x, y, year) %>%
  summarise(AOT40 = sum(AOT40) * 1e-3, .groups = "drop") %>%
  inner_join(fraction) %>%
  as_tibble()

ozone_ctr <- read_rds("data/ozone_cft.rds") %>%
  filter(peak_level == 60) %>%
  pull() %>%
  .[[1]]

ozone_ctr <- cldr %>%
  lazy_dt() %>%
  inner_join(ozone_ctr, by = c("x", "y", "month")) %>%
  group_by(crop, year, x, y) %>%
  summarise(AOT40 = sum(AOT40_cft) * 1e-3, .groups = "drop") %>%
  as_tibble()


data <-
  read_xlsx("data/inputs/Field_studies/FengNF2022.xlsx",
    sheet = "AOT40-RY", range = "A4:M1000"
  ) %>%
  janitor::remove_empty(c("rows", "cols")) %>%
  rename(AOT40 = `AOT40 (ppm h)`)

models <-
  list(
    c("RY_maize", "RY_maize_LB95", "RY_maize_UB95"),
    c("RY_inbred_rice", "RY_inbred_rice_LB95", "RY_inbred_rice_UB95"),
    c("RY_inbred_rice", "RY_inbred_rice_LB95", "RY_inbred_rice_UB95"),
    c("RY_wheat", "RY_wheat_LB95", "RY_wheat_UB95")
  ) %>%
  map(function(acrops) {
    acrops %>%
      str_c(" ~ s(AOT40)") %>%
      map(as.formula) %>%
      map(~ gam(.x, data = data)) %>%
      `names<-`(c("50%", "95%", "5%"))
  })

ozone <- ozone %>% nest(fdata = -crop)
ozone_ctr <- ozone_ctr %>% nest(fdata = -crop)

results <-
  tibble(
    crop = ozone$crop,
    results = pro_pmap(
      list(
        ozone$fdata, ozone_ctr$fdata, models
      ),
      function(aozone, aozone_ctr, amodel) {
        oh <- aozone %>%
          bind_cols(map_dfc(amodel, predict, newdata = ., cluster = cl)) %>%
          select(-AOT40) %>%
          rename_with(~ str_c("hist", .x), contains("%"))

        oc <- aozone_ctr %>%
          bind_cols(map_dfc(amodel, predict, newdata = ., cluster = cl)) %>%
          select(-AOT40) %>%
          rename_with(~ str_c("ctr", .x), contains("%"))

        inner_join(oh, oc) %>%
          mutate(
            `50%` = `ctr50%` - `hist50%`,
            `5%` = `ctr5%` - `hist5%`,
            `95%` = `ctr95%` - `hist95%`,
            .keep = "unused"
          )
      }
    )
  )

qsave(results, "data/impacts_ozone_Feng2022.qs", nthreads = qn)
