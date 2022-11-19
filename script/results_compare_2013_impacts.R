source("script/loadPackages.R")
source("script/loadFunctions.R")

cldr <- read_rds("data/outputs/calendar/tidied.rds") %>%
  filter((MA - `GR&EM`) >= 2) %>%
  trim_xy()

fraction <-
  tibble(
    crop = c("Maize", "Rice(LR)", "Rice(SR&ER)", "Wheat"),
    data = map(c("Maize", "Rice", "Rice", "Wheat"), function(acrop) {
      afile <- str_c("data/outputs/masks/mask_", acrop, ".tif")
      rast(afile) %>%
        as.data.frame(xy = T, na.rm = F) %>%
        pivot_longer(
          -c(x, y),
          names_to = "year",
          names_transform = list(year = as.integer),
          values_drop_na = T,
          values_to = "fraction"
        )
    })
  ) %>%
  unnest() %>%
  trim_xy()

f1 <- read_rds("data/boots_f1.rds") %>%
  mutate(crop = crop_parent) %>%
  bind_rows((.) %>% filter(crop == "Rice") %>% mutate(crop = "Rice(LR)")) %>%
  mutate(crop = fifelse(crop == "Rice", "Rice(SR&ER)", crop))


# AOD ---------------------------------------------------------------------

AOD <- read_rds("data/outputs/aerosol/AOD.rds")

aod_ctr <- map_dfr(unique(AOD$year), function(ayear) {
  AOD %>%
    filter(year == 2013) %>%
    mutate(year = ayear)
})

AOD <- cldr %>%
  lazy_dt() %>%
  inner_join(AOD) %>%
  group_by(crop, x, y, year) %>%
  summarise(AOD = mean(AOD), .groups = "drop") %>%
  inner_join(fraction) %>%
  as_tibble() %>%
  nest(fdata = -c(crop, year)) %>%
  arrange(crop, year)

aod_ctr <- cldr %>%
  lazy_dt() %>%
  inner_join(aod_ctr) %>%
  group_by(crop, x, y, year) %>%
  summarise(AOD_cft = mean(AOD), .groups = "drop") %>%
  as_tibble() %>%
  nest(aod_data = -c(crop, year)) %>%
  arrange(crop, year)

impacts <-
  reduce(list(AOD, f1, aod_ctr), inner_join)

# adata = impacts$fdata[[50]]
# coefs = impacts$coefs[[50]]
# aod_data = impacts$aod_data[[50]]

plan(multisession, workers = 5)

impacts <- impacts %>%
  mutate(AOD_results = future_pmap(list(fdata, coefs, aod_data), function(adata, coefs, aod_data) {
    coefs <- coefs %>%
      select(id, term, estimate) %>%
      pivot_wider(names_from = term, values_from = estimate) %>%
      select(contains("AOD"))

    cal_data <- adata %>%
      inner_join(aod_data, by = c("x", "y"))

    rel_X <- cal_data %>%
      mutate(AOD1 = AOD_cft - AOD, AOD2 = AOD_cft^2 - AOD^2) %>%
      select(AOD1, AOD2)

    # cal relative impact
    rel_results <- tcrossprod(as.matrix(rel_X), as.matrix(coefs)) %>%
      expm1() %>%
      rowQuantiles(probs = c(0.05, 0.5, 0.95)) %>%
      as_tibble() %>%
      bind_cols(cal_data %>% select(x, y, fraction))

    return(rel_results)
  }, .progress = T)) %>%
  select(-c(fdata, coefs, aod_data))

plan(sequential)

qsave(impacts, "data/impacts_2013_aerosol.qs", nthread = qn)


# ozone ---------------------------------------------------------------------

cldr <-
  bind_rows(
    cldr %>%
      filter(crop == "Maize"),
    cldr %>%
      filter(str_detect(crop, "Rice") & month - `GR&EM` != 6),
    cldr %>%
      filter(crop == "Wheat" & MA != month)
  )


ozone <- read_rds("data/outputs/ozone/tidied.rds")

ozone_ctr <- map_dfr(unique(ozone$year), function(ayear) {
  ozone %>%
    filter(year == 2013) %>%
    mutate(year = ayear)
})

ozone <- cldr %>%
  lazy_dt() %>%
  inner_join(ozone) %>%
  group_by(crop, x, y, year) %>%
  summarise(AOT40 = sum(AOT40), .groups = "drop") %>%
  inner_join(fraction) %>%
  as_tibble() %>%
  nest(fdata = -c(crop, year)) %>%
  arrange(crop, year)

ozone_ctr <- cldr %>%
  lazy_dt() %>%
  inner_join(ozone_ctr) %>%
  group_by(crop, x, y, year) %>%
  summarise(AOT40_cft = sum(AOT40), .groups = "drop") %>%
  as_tibble() %>%
  nest(ozone_data = -c(crop, year)) %>%
  arrange(crop, year)

impacts <-
  reduce(list(ozone, f1, ozone_ctr), inner_join)

# adata = impacts$fdata[[9]]
# coefs = impacts$coefs[[9]]
# ozone_data = impacts$ozone_data[[9]]

plan(multisession, workers = 5)

impacts <- impacts %>%
  mutate(ozone_results = future_pmap(list(fdata, coefs, ozone_data), function(adata, coefs, ozone_data) {
    coefs <- coefs %>%
      select(id, term, estimate) %>%
      pivot_wider(names_from = term, values_from = estimate) %>%
      select(contains("AOT40"))

    cal_data <- adata %>%
      inner_join(ozone_data, by = c("x", "y"))

    rel_X <- cal_data %>%
      mutate(AOT401 = AOT40_cft - AOT40) %>%
      select(AOT401)

    # cal relative impact
    rel_results <- tcrossprod(as.matrix(rel_X), as.matrix(coefs)) %>%
      expm1() %>%
      rowQuantiles(probs = c(0.05, 0.5, 0.95)) %>%
      as_tibble() %>%
      bind_cols(cal_data %>% select(x, y, fraction))

    return(rel_results)
  }, .progress = T)) %>%
  select(-c(fdata, coefs, ozone_data))

plan(sequential)

qsave(impacts, "data/impacts_2013_AOT40.qs", nthread = qn)

# plot --------------------------------------------------------------------

order <- c(
  "North China",
  "Northeast China",
  "East China",
  "China",
  "South China",
  "Central China",
  "Southwest China",
  "Northwest China"
)

province_plot <- st_read("https://geo.datav.aliyun.com/areas_v3/bound/100000_full.json") %>%
  mutate(
    region = fcase(
      name %in% c("辽宁省", "吉林省", "黑龙江省"),
      "Northeast China",
      name %in% c("上海市", "江苏省", "浙江省", "安徽省", "福建省", "江西省", "山东省", "台湾省"),
      "East China",
      name %in% c("北京市", "天津市", "河北省", "山西省", "内蒙古自治区"),
      "North China",
      name %in% c("河南省", "湖北省", "湖南省"),
      "Central China",
      name %in% c("广东省", "广西壮族自治区", "海南省", "香港特别行政区", "澳门特别行政区"),
      "South China",
      name %in% c("四川省", "贵州省", "云南省", "西藏自治区", "重庆市"),
      "Southwest China",
      name %in% c("陕西省", "甘肃省", "青海省", "宁夏回族自治区", "新疆维吾尔自治区"),
      "Northwest China"
    )
  )

region_plot <-
  province_plot %>%
  st_make_valid() %>%
  group_by(region) %>%
  summarise(geometry = st_union(geometry))

plot(region_plot)

province <-
  st_read("data/inputs/shp/2020省矢量.shp") %>%
  select(-year) %>%
  mutate(
    region = fcase(
      省 %in% c("辽宁省", "吉林省", "黑龙江省"),
      "Northeast China",
      省 %in% c("上海市", "江苏省", "浙江省", "安徽省", "福建省", "江西省", "山东省", "台湾省"),
      "East China",
      省 %in% c("北京市", "天津市", "河北省", "山西省", "内蒙古自治区"),
      "North China",
      省 %in% c("河南省", "湖北省", "湖南省"),
      "Central China",
      省 %in% c("广东省", "广西壮族自治区", "海南省", "香港特别行政区", "澳门特别行政区"),
      "South China",
      省 %in% c("四川省", "贵州省", "云南省", "西藏自治区", "重庆市"),
      "Southwest China",
      省 %in% c("陕西省", "甘肃省", "青海省", "宁夏回族自治区", "新疆维吾尔自治区"),
      "Northwest China"
    )
  )

region <-
  bind_rows(
    province %>%
      group_by(region) %>%
      summarise(geometry = st_union(geometry)),
    province %>%
      summarise(geometry = st_union(geometry)) %>%
      mutate(region = "China")
  )

design <- "
H#A##
#####
####B
GF##C
##E#D
"

aerosol <- qread("data/impacts_2013_aerosol.qs", nthreads = qn)

aerosol <-
  aerosol %>%
  mutate(AOD_results = pro_map(AOD_results, function(adata) {
    temp <- adata %>%
      relocate(x, y) %>%
      rast(type = "xyz", crs = "epsg:4326")

    fraction <- temp[["fraction"]]
    temp <- temp["%"]

    province_level <-
      exact_extract(temp,
        province,
        "weighted_mean",
        weights = fraction,
        progress = F
      ) %>%
      bind_cols(province, .) %>%
      st_drop_geometry() %>%
      pivot_longer(contains("%"),
        names_prefix = "weighted_mean.",
        names_to = "name"
      ) %>%
      pivot_wider()

    region_level <- exact_extract(temp,
      region,
      "weighted_mean",
      weights = fraction,
      progress = F
    ) %>%
      bind_cols(region, .) %>%
      st_drop_geometry() %>%
      pivot_longer(contains("%"),
        names_prefix = "weighted_mean.",
        names_to = "name"
      ) %>%
      pivot_wider()

    lst(province_level, region_level)
  }))

ozone <- qread("data/impacts_2013_AOT40.qs", nthreads = qn)

ozone <-
  ozone %>%
  mutate(ozone_results = pro_map(ozone_results, function(adata) {
    temp <- adata %>%
      relocate(x, y) %>%
      rast(type = "xyz", crs = "epsg:4326")

    fraction <- temp[["fraction"]]
    temp <- temp["%"]

    province_level <-
      exact_extract(temp,
        province,
        "weighted_mean",
        weights = fraction,
        progress = F
      ) %>%
      bind_cols(province, .) %>%
      st_drop_geometry() %>%
      pivot_longer(contains("%"),
        names_prefix = "weighted_mean.",
        names_to = "name"
      ) %>%
      pivot_wider()

    region_level <- exact_extract(temp,
      region,
      "weighted_mean",
      weights = fraction,
      progress = F
    ) %>%
      bind_cols(region, .) %>%
      st_drop_geometry() %>%
      pivot_longer(contains("%"),
        names_prefix = "weighted_mean.",
        names_to = "name"
      ) %>%
      pivot_wider()

    lst(province_level, region_level)
  }))

p4 <- aerosol %>%
  unnest_wider(AOD_results) %>%
  select(-province_level) %>%
  unnest() %>%
  group_by(crop_parent, region, year) %>%
  summarise(across(contains("%"), fmean), .groups = "drop") %>%
  mutate(
    region = factor(region, order),
    across(contains("%"), ~ .x * 1e2)
  ) %>%
  ggplot(aes(
    x = year,
    y = `50%`,
    color = factor(region, order),
    alpha = factor(region, order)
  ))

p5 <- ozone %>%
  unnest_wider(ozone_results) %>%
  select(-province_level) %>%
  unnest() %>%
  group_by(crop_parent, region, year) %>%
  summarise(across(contains("%"), fmean), .groups = "drop") %>%
  mutate(
    region = factor(region, order),
    across(contains("%"), ~ .x * 1e2)
  ) %>%
  ggplot(aes(
    x = year,
    y = `50%`,
    color = factor(region, order),
    alpha = factor(region, order)
  ))


patch <- wrap_plots(p4, p5, guides = "collect") +
  plot_annotation(tag_levels = "a") &
  facet_wrap(~crop_parent, ncol = 1, scales = "free_y") &
  geom_hline(yintercept = 0, linetype = "dashed") &
  geom_line() &
  geom_point(size = 2) &
  stat_poly_eq(use_label(c("eq.label", "p.value.label")),
    data = . %>% filter(region == "China"),
    family = "Roboto Condensed"
  ) &
  theme_half_open(16, font_family = "Roboto Condensed") &
  background_grid() &
  scale_color_manual(
    name = NULL,
    values = MetBrewer::met.brewer("Juarez", n = 7) %>% append("black", 3)
  ) &
  scale_alpha_manual(name = NULL, values = rep(0.2, 7) %>% append(1, 3)) &
  xlab("Year") &
  ylab("Percentage change in SIF")

ggsave(
  "figures/hist_impact_2013.pdf",
  patch,
  width = 4,
  height = 3,
  scale = 3
)
