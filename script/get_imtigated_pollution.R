source("script/loadPackages.R")
source("script/loadFunctions.R")

# aerosol -----------------------------------------------------------------

pm <- read_rds("data/outputs/aerosol/PM25.rds")
aod <- read_rds("data/outputs/aerosol/AOD.rds")

joined <- inner_join(pm, aod) %>%
  fgroup_by(x, y, year) %>%
  fsummarise(across(c(PM25, AOD), fmean))

model <- feols(AOD ~ PM25 | x^y,
  data = joined,
  cluster = ~ x^y, nthreads = 0, combine.quick = F
)

esttable(model)

month_climatology <- aod %>%
  fgroup_by(x, y, month) %>%
  fsummarise(AOD_month = fmean(AOD)) %>%
  fgroup_by(x, y) %>%
  fmutate(AOD_month_percent = AOD_month / fmean(AOD_month), .keep = "unused") %>%
  fungroup()

aod_cft_all <-
  tibble(
    pm_level = seq(0, 90, by = 5),
    aod_data = map(pm_level, function(anum) {
      month_climatology %>%
        distinct(x, y) %>%
        mutate(PM25 = anum) %>%
        mutate(
          AOD = predict(model, .),
          AOD = fifelse(AOD < 0, 0, AOD)
        ) %>%
        select(-PM25) %>%
        inner_join(month_climatology) %>%
        mutate(AOD_cft = AOD_month_percent * AOD, .keep = "unused")
    })
  )

saveRDS(aod_cft_all, "data/aod_cft.rds")

# aerosol2 -----------------------------------------------------------------

pm <- read_rds("data/outputs/aerosol/PM25.rds")
aod <- read_rds("data/outputs/aerosol/AOD.rds")

joined <- inner_join(pm, aod)

model <- feols(AOD ~ PM25 | x^y + month,
  data = joined,
  cluster = ~ x^y, nthreads = 0, combine.quick = F
)

esttable(model)

month_climatology <- pm %>%
  fgroup_by(x, y, month) %>%
  fsummarise(PM25_month = fmean(PM25)) %>%
  fgroup_by(x, y) %>%
  fmutate(PM25_month_percent = PM25_month / fmean(PM25_month), .keep = "unused") %>%
  fungroup()

aod_cft_all <-
  tibble(
    pm_level = seq(0, 50, by = 5),
    aod_data = map(pm_level, function(anum) {
      month_climatology %>%
        mutate(PM25 = anum * PM25_month_percent) %>%
        mutate(
          AOD_cft = predict(model, .),
          AOD_cft = fifelse(AOD_cft < 0, 0, AOD_cft)
        ) %>%
        select(-starts_with("PM"))
    })
  )

saveRDS(aod_cft_all, "data/aod_cft2.rds")


# ozone -------------------------------------------------------------------

ozone <- read_rds("data/outputs/ozone/tidied.rds") # unit is ug/m3

yearly <- ozone %>%
  lazy_dt() %>%
  group_by(x, y, year) %>%
  summarise(peak = fmean(O3[month %in% 4:9]), annual = fmean(O3), .groups = "drop") %>%
  as_tibble()

model <- feols(annual ~ peak | x^y, yearly, cluster = ~ x^y, nthreads = 0, combine.quick = F)

month_climatology <- ozone %>%
  fgroup_by(x, y, month) %>%
  fsummarise(O3_month = fmean(O3)) %>%
  fgroup_by(x, y) %>%
  fmutate(O3_month_percent = O3_month / fmean(O3_month), .keep = "unused") %>%
  fungroup()

veg_model <- read_rds("data/outputs/ozone/gam.rds")

require(parallel)
cl <- makeCluster(qn)
ozone_cft_all <-
  tibble(
    peak_level = seq(30, 150, 5),
    ozone_data = pro_map(peak_level, function(anum) {
      month_climatology %>%
        distinct(x, y) %>%
        mutate(peak = anum) %>%
        mutate(
          O3 = predict(model, .),
          O3 = fifelse(O3 < 0, 0, O3)
        ) %>%
        select(-peak) %>%
        inner_join(month_climatology) %>%
        mutate(O3 = O3_month_percent * O3, .keep = "unused") %>%
        fmutate(
          W126_cft = exp(predict(veg_model$model_w126, ., cluster = cl)),
          AOT40_cft = exp(predict(veg_model$model_aot40, ., cluster = cl))
        ) %>%
        select(-O3)
    })
  )

saveRDS(ozone_cft_all, "data/ozone_cft.rds")
