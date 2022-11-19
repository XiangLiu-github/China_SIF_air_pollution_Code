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
  as_tibble() %>%
  nest(fdata = -c(crop, year)) %>%
  arrange(crop, year)

ozone_ctr <- read_rds("data/ozone_cft.rds")

ozone_ctr <- ozone_ctr %>%
  mutate(ozone_data = pro_map(ozone_data, function(adata) {
    cldr %>%
      lazy_dt() %>%
      inner_join(adata, by = c("x", "y", "month")) %>%
      group_by(crop, year, x, y) %>%
      summarise(AOT40_cft = sum(AOT40_cft) * 1e-3, .groups = "drop") %>%
      as_tibble()
  })) %>%
  unnest() %>%
  nest(ozone_data = -c(crop, year))

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

models = tibble(crop = c('Maize', 'Rice(LR)', 'Rice(SR&ER)', 'Wheat'),
                model = models)

impacts <-
  reduce(list(ozone, models, ozone_ctr), inner_join)

# adata = impacts$fdata[[1]]
# coefs = impacts$model[[1]]
# ozone_data = impacts$ozone_data[[1]]

impacts <- impacts %>%
  mutate(ozone_results = pro_pmap(list(fdata, model, ozone_data), function(adata, coefs, ozone_data) {
    
    oh <- adata %>%
      bind_cols(map_dfc(coefs, predict, newdata = ., cluster = cl)) %>%
      select(-AOT40) %>%
      rename_with(~ str_c("hist", .x), contains("%"))
    
    oc <- ozone_data %>%
      rename(AOT40 = AOT40_cft) %>% 
      bind_cols(map_dfc(coefs, predict, newdata = ., cluster = cl)) %>%
      select(-AOT40) %>%
      rename_with(~ str_c("ctr", .x), contains("%"))
    
    rel_results = inner_join(oh, oc) %>%
      mutate(
        `50%` = `ctr50%` - `hist50%`,
        `5%` = `ctr5%` - `hist5%`,
        `95%` = `ctr95%` - `hist95%`,
        .keep = "unused"
      )
    
    return(rel_results)
  })) %>%
  select(-c(fdata, model, ozone_data))

qsave(impacts, "data/impacts_ozone_Feng2022.qs", nthreads = qn)


















