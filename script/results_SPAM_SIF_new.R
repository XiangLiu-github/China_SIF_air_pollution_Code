source("script/loadPackages.R")
source("script/loadFunctions.R")

gs <- read_rds("data/outputs/calendar/tidied.rds") %>%
  trim_xy()

fraction <- map(c("Maize", "Rice", "Wheat"), ~ rast(str_c("data/outputs/masks/mask_", .x, ".tif"))) %>%
  rast() %>%
  rep(each = 3)

names(fraction) <- 1:nlyr(fraction)

GOSIF <- read_rds("data/outputs/SIF/GOSIF.rds")

GOSIF <- inner_join(gs, GOSIF) %>%
  fgroup_by(crop, x, y, year) %>%
  fsummarise(GOSIF_mean = mean(GOSIF), GOSIF_peak = max(GOSIF), GOSIF_sum = sum(GOSIF))

GOSIF <- GOSIF %>%
  mutate(crop = fifelse(str_detect(crop, "Rice"), "Rice", crop)) %>%
  fgroup_by(crop, x, y, year) %>%
  fsummarise(across(c(GOSIF_mean, GOSIF_peak, GOSIF_sum), fmean))

# GOSIF %>%
#   filter(crop == 'Rice', year == 2015) %>%
#   select(x, y, starts_with('GOSIF')) %>%
#   rast() %>%
#   plot()

GOSIF <- GOSIF %>%
  pivot_longer(starts_with("GOSIF")) %>%
  pivot_wider(names_from = c(crop, year, name), names_sep = "_") %>%
  rast(crs = "epsg:4326")

GOSIF <- GOSIF[[sort(names(GOSIF))]]

NIRv <- read_rds("data/outputs/NIRv/tidied.rds") %>%
  select(x, y, year, month, NIRvc)

NIRv <- inner_join(gs, NIRv) %>%
  fgroup_by(crop, x, y, year) %>%
  fsummarise(NIRv_mean = mean(NIRvc), NIRv_peak = max(NIRvc), NIRv_sum = sum(NIRvc))

NIRv <- NIRv %>%
  mutate(crop = fifelse(str_detect(crop, "Rice"), "Rice", crop)) %>%
  fgroup_by(crop, x, y, year) %>%
  fsummarise(across(c(NIRv_mean, NIRv_peak, NIRv_sum), fmean))

NIRv <- NIRv %>%
  pivot_longer(starts_with("NIRv")) %>%
  pivot_wider(names_from = c(crop, year, name), names_sep = "_") %>%
  rast(crs = "epsg:4326")

NIRv <- NIRv[[sort(names(NIRv))]]


# SPAM --------------------------------------------------------------------

shp <- st_read("../tidy_crop_statistics/data/inputs/shp/2020年县.shp") %>% select(-year)

GOSIF_extracted <- exact_extract(GOSIF, shp, "weighted_mean", weights = fraction) %>%
  bind_cols(shp %>% st_drop_geometry())

GOSIF_extracted <- GOSIF_extracted %>%
  pivot_longer(starts_with("weighted_mean."),
    names_prefix = "weighted_mean.", names_sep = "_", names_to = c("crop", "year", "name", "summary"), names_transform = list(year = as.integer)
  ) %>%
  pivot_wider(values_fn = sum) %>%
  separate(summary, c("summary", NA), "\\.")

spam_all <- map2_dfr(
  c("Wheat", "Rice", "Maize"),
  c("WHEA", "RICE", "MAIZ"),
  function(acrop1, acrop) {
    temp <- rast(str_c("/vsizip/../data_archive/SPAM/spam2000v3.0.7_global_harvested-area.geotiff.zip/spam2000V3r107_global_H_", acrop, "_A.tif"))

    harv <- c(
      str_c("/vsizip/../data_archive/SPAM/spam2000v3.0.7_global_harvested-area.geotiff.zip/spam2000V3r107_global_H_", acrop, "_A.tif"),
      str_c("/vsizip/../data_archive/SPAM/spam2005v3r2_global_harv_area.geotiff.zip/geotiff_global_harv_area/SPAM2005V3r2_global_H_TA_", acrop, "_A.tif"),
      str_c("/vsizip/../data_archive/SPAM/spam2010v2r0_global_harv_area.geotiff.zip/spam2010V2r0_global_H_", acrop, "_A.tif")
    ) %>%
      map(rast) %>%
      map(~ extend(.x, temp)) %>%
      rast() %>%
      `names<-`(str_c("harv_", seq(2000, 2010, by = 5)))

    prod <- c(
      str_c("/vsizip/../data_archive/SPAM/spam2000v3.0.7_global_production.geotiff.zip/spam2000V3r107_global_R_", acrop, "_A.tif"),
      str_c("/vsizip/../data_archive/SPAM/spam2005v3r2_global_prod.geotiff.zip/geotiff_global_prod/SPAM2005V3r2_global_P_TA_", acrop, "_A.tif"),
      str_c("/vsizip/../data_archive/SPAM/spam2010v2r0_global_prod.geotiff.zip/spam2010V2r0_global_P_", acrop, "_A.tif")
    ) %>%
      map(rast) %>%
      map(~ extend(.x, temp)) %>%
      rast() %>%
      `names<-`(str_c("prod_", seq(2000, 2010, by = 5)))


    exact_extract(c(harv, prod), shp, "sum") %>%
      bind_cols(shp %>% st_drop_geometry()) %>%
      pivot_longer(starts_with("sum"),
        names_prefix = "sum.", names_to = c("name", "year"), names_sep = "_", names_transform = list(year = as.integer)
      ) %>%
      pivot_wider(values_fn = sum) %>%
      mutate(crop = acrop1)
  }
)

fraction_extracted <- map(c("Maize", "Rice", "Wheat"), ~ rast(str_c("data/outputs/masks/mask_", .x, ".tif"))) %>%
  rast() %>%
  `names<-`(str_c(rep(c("Maize", "Rice", "Wheat"), each = 20), "_", names(.))) %>%
  exact_extract(shp, "sum") %>%
  bind_cols(shp %>% st_drop_geometry()) %>%
  pivot_longer(starts_with("sum."),
    names_to = c("crop", "year"), names_sep = "_", names_transform = list(year = as.integer), names_prefix = "sum.",
    values_to = "frac"
  )

spam_all <- spam_all %>%
  filter(!(crop %in% c("Wheat", "Rice") & year == 2000))

list(GOSIF_extracted, fraction_extracted, spam_all) %>%
  reduce(inner_join) %>%
  slice_max(frac, prop = 0.6, by = crop) %>%
  filter(省级 == "河南省") %>%
  ggplot(aes(x = GOSIF, y = prod / harv, color = factor(year))) +
  # ggplot(aes(x = frac, y = harv, color = factor(year))) +
  facet_grid2(vars(crop), vars(summary), scales = "free", independent = "all") +
  geom_point() +
  geom_smooth(method = "lm") +
  stat_poly_eq()

diffs <- list(GOSIF_extracted, fraction_extracted, spam_all) %>%
  reduce(inner_join) %>%
  group_by(
    地名, 区划码, 县级, 县级码, 县级类, 地级, 地级码,
    地级类, 省级, 省级码, 省级类, 曾用名, 备注, ENG_NAME, VAR_NAME, code, NAME_3,
    VAR_NAME3, GID_3, TYPE_3, NAME_2, VAR_NAME2, GID_2, TYPE_2, NAME_1,
    VAR_NAME1, GID_1, TYPE_1, crop, summary
  ) %>%
  relocate(year) %>%
  arrange(year, .by_group = T) %>%
  mutate(
    yield = prod / harv,
    GOSIF = log(GOSIF) - lag(log(GOSIF), order_by = year),
    yield = log(yield) - lag(log(yield), order_by = year)
  ) %>%
  ungroup()

diffs %>%
  slice_max(frac, prop = 0.6, by = crop) %>%
  # filter(省级 == '河南省') %>%
  ggplot(aes(x = GOSIF, y = yield, color = factor(year))) +
  # ggplot(aes(x = frac, y = harv, color = factor(year))) +
  facet_grid2(vars(crop), vars(summary), scales = "free", independent = "all") +
  geom_point() +
  geom_smooth(method = "lm") +
  stat_poly_eq()

# province level ----------------------------------------------------------

shp <- st_read("data/inputs/shp/2020省矢量.shp") %>% select(-year)

GOSIF_extracted <- exact_extract(GOSIF, shp, "weighted_mean", weights = fraction) %>%
  bind_cols(shp %>% st_drop_geometry()) %>%
  pivot_longer(starts_with("weighted_mean."),
    names_prefix = "weighted_mean.", names_sep = "_", names_to = c("crop", "year", "name", "summary"), names_transform = list(year = as.integer)
  ) %>%
  pivot_wider(values_fn = sum) %>%
  separate(summary, c("summary", NA), "\\.")

NIRv_extracted <- exact_extract(NIRv, shp, "weighted_mean", weights = fraction) %>%
  bind_cols(shp %>% st_drop_geometry()) %>%
  pivot_longer(starts_with("weighted_mean."),
    names_prefix = "weighted_mean.", names_sep = "_", names_to = c("crop", "year", "name", "summary"), names_transform = list(year = as.integer)
  ) %>%
  pivot_wider(values_fn = sum) %>%
  separate(summary, c("summary", NA), "\\.")

fraction_extracted <- map(c("Maize", "Rice", "Wheat"), ~ rast(str_c("data/outputs/masks/mask_", .x, ".tif"))) %>%
  rast() %>%
  `names<-`(str_c(rep(c("Maize", "Rice", "Wheat"), each = 20), "_", names(.))) %>%
  exact_extract(shp, "sum") %>%
  bind_cols(shp %>% st_drop_geometry()) %>%
  pivot_longer(starts_with("sum."),
    names_to = c("crop", "year"), names_sep = "_", names_transform = list(year = as.integer), names_prefix = "sum.",
    values_to = "frac"
  )

Area <-
  tibble(
    files = list.files("data/inputs/province_data/area/", full.names = T),
    data = map(files, ~ fread(.x, skip = 3, nrows = 31)),
    crop = map_chr(files, ~ fread(.x, skip = 1, nrows = 1) %>% pull(V1))
  ) %>%
  select(-files) %>%
  mutate(crop = str_remove(crop, "指标：") %>% str_remove("播种面积\\(千公顷\\)")) %>%
  unnest() %>%
  rename(year = 时间) %>%
  pivot_longer(-c(crop, year), names_to = "province", values_to = "area") %>%
  drop_na() %>%
  mutate(year = parse_number(year)) %>%
  mutate(area_mean = mean(area), .by = c(crop, province))

Production <-
  tibble(
    files = list.files("data/inputs/province_data/production/", full.names = T),
    data = map(files, ~ fread(.x, skip = 3, nrows = 31)),
    crop = map_chr(files, ~ fread(.x, skip = 1, nrows = 1) %>% pull(V1))
  ) %>%
  select(-files) %>%
  mutate(crop = str_remove(crop, "指标：") %>% str_remove("产量\\(万吨\\)")) %>%
  unnest() %>%
  rename(year = 时间) %>%
  pivot_longer(-c(crop, year), names_to = "province", values_to = "production") %>%
  drop_na() %>%
  mutate(year = parse_number(year))

stats <-
  inner_join(Area, Production) %>%
  mutate(region = fcase(
    province %in% c("辽宁省", "吉林省", "黑龙江省"), "Northeast China",
    province %in% c("上海市", "江苏省", "浙江省", "安徽省", "福建省", "江西省", "山东省", "台湾省"), "East China",
    province %in% c("北京市", "天津市", "河北省", "山西省", "内蒙古自治区"), "North China",
    province %in% c("河南省", "湖北省", "湖南省"), "Central China",
    province %in% c("广东省", "广西壮族自治区", "海南省", "香港特别行政区", "澳门特别行政区"), "South China",
    province %in% c("四川省", "贵州省", "云南省", "西藏自治区", "重庆市"), "Southwest China",
    province %in% c("陕西省", "甘肃省", "青海省", "宁夏回族自治区", "新疆维吾尔自治区"), "Northwest China"
  )) %>%
  mutate(production = production * 1e4) %>%
  mutate(
    yield = production / area,
    crop = fcase(
      crop == "小麦", "Wheat",
      crop == "玉米", "Maize",
      crop == "稻谷", "Rice"
    )
  ) %>%
  drop_na() %>%
  filter(year %in% 2001:2019) %>%
  rename(省 = province)

# anti_join(stats, GOSIF_extracted)

# list(GOSIF_extracted, fraction_extracted, stats) %>%
#   reduce(inner_join) %>%
#   mutate(across(c(GOSIF, yield), fscale), .by = c(省, crop, summary)) %>%
#   summarise(cor = weights::wtd.cors(yield, GOSIF, area_mean)[1,1], .by = c(crop, summary))

p1 <- list(GOSIF_extracted, fraction_extracted, stats) %>%
  reduce(inner_join) %>%
  mutate(across(c(GOSIF, yield), fscale), .by = c(省, crop, summary)) %>%
  arrange(area_mean) %>%
  filter(area_mean >= 10) %>%
  ggplot(aes(x = GOSIF, y = yield)) +
  facet_grid2(vars(str_c("GOSIF ", summary)), vars(crop)) +
  geom_abline(slope = 1, intercept = 0, color = "blue") +
  geom_pointdensity(aes(color = stat(ndensity)), size = 2) +
  geom_smooth(
    method = "lm", mapping = aes(weight = log(area_mean)),
    color = "red", show.legend = FALSE, se = F
  ) +
  geom_text_repel(aes(x = -3, y = 3, label = sprintf("%0.2f", round(cor, digits = 2))), inherit.aes = F, data = . %>% summarise(cor = weights::wtd.cors(yield, GOSIF, log(area_mean))[1, 1], .by = c(crop, summary)), family = "Roboto Condensed", size = 5) +
  scale_color_gradientn(colours = PNWColors::pnw_palette("Moth") %>% rev()) +
  ylab("Normalized yield") +
  xlab("Normalized SIF") +
  theme_half_open(18, font_family = "Roboto Condensed") +
  background_grid()

ggsave("figures/province_SIF_validation.pdf", p1, width = 1, height = 1, scale = 10)

p2 <- list(NIRv_extracted, fraction_extracted, stats) %>%
  reduce(inner_join) %>%
  mutate(across(c(NIRv, yield), fscale), .by = c(省, crop, summary)) %>%
  arrange(area_mean) %>%
  filter(area_mean >= 10) %>%
  ggplot(aes(x = NIRv, y = yield)) +
  facet_grid2(vars(str_c("NIRv ", summary)), vars(crop)) +
  geom_abline(slope = 1, intercept = 0, color = "blue") +
  geom_pointdensity(aes(color = stat(ndensity)), size = 2) +
  geom_smooth(
    method = "lm", mapping = aes(weight = log(area_mean)),
    color = "red", show.legend = FALSE, se = F
  ) +
  geom_text_repel(aes(x = -3, y = 3, label = sprintf("%0.2f", round(cor, digits = 2))), inherit.aes = F, data = . %>% summarise(cor = weights::wtd.cors(yield, NIRv, log(area_mean))[1, 1], .by = c(crop, summary)), family = "Roboto Condensed", size = 5) +
  scale_color_gradientn(colours = PNWColors::pnw_palette("Moth") %>% rev()) +
  ylab("Normalized yield") +
  xlab("Normalized NIRv") +
  theme_half_open(18, font_family = "Roboto Condensed") +
  background_grid()

ggsave("figures/province_NIRv_validation.pdf", p2, width = 1, height = 1, scale = 10)

detrend = function(value, year){

  lm(value ~ year, data = NULL)$residuals
  
}

list(NIRv_extracted, fraction_extracted, stats) %>%
  reduce(inner_join) %>%
  drop_na() %>% 
  # mutate(across(c(NIRv, yield), ~ detrend(.x, year)), .by = c(省, crop, summary)) %>%
  arrange(area_mean) %>%
  filter(area_mean >= 10) %>%
  summarise(cor = cor(NIRv, yield), 
            cor_f = cor(frac, area),
            .by = c(ENG_NAME, crop, summary, area_mean)) %>% 
  # filter(cor_f >= 0.5) %>% 
  ggplot(aes(x = area_mean, y = cor)) +
  facet_grid(vars(summary), vars(crop)) +
  geom_point() +
  geom_text_repel(aes(label = ENG_NAME))

#
# fpro = stats %>%
#   distinct(crop, 省, area_mean) %>%
#   slice_max(area_mean, n = 9, by = crop)
#
# p1 = list(GOSIF_extracted, fraction_extracted, stats) %>%
#   reduce(inner_join) %>%
#   mutate(across(c(GOSIF, yield), fscale), .by = c(省, crop, summary)) %>%
#   arrange(area_mean) %>%
#   inner_join(fpro) %>%
#   filter(summary == 'mean', crop == 'Maize') %>%
#   mutate(ENG_NAME = fct_reorder(ENG_NAME, area_mean, .desc = T)) %>%
#   ggplot(aes(x = year)) +
#   facet_wrap(~ ENG_NAME, scales = "free_y") +
#   geom_line(aes(y = GOSIF), color = 'red') +
#   geom_point(aes(y = GOSIF), color = 'red', size = 2) +
#   geom_line(aes(y = yield)) +
#   geom_point(aes(y = yield)) +
#   geom_text(aes(x = 2002, y = 2, label = sprintf("%0.2f", round(cor, digits = 2))), inherit.aes = F, data = . %>% summarise(cor = cor(GOSIF, yield, use = "complete.obs"), .by = c(crop, ENG_NAME, summary)), check_overlap = T, family = "Roboto Condensed", size = 5) +
#   ylab('Normalized SIF or yield') +
#   xlab('Year') +
#   theme_half_open(18, font_family = "Roboto Condensed") +
#   background_grid()
#
# ggsave('test.pdf', p1, width = 1.5, height = 1, scale = 8)
#
