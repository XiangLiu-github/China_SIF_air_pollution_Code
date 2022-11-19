source("script/loadPackages.R")
source("script/loadFunctions.R")

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
  mutate(year = parse_number(year))

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
  mutate(production = production * 1e4)

stats_region <-
  stats %>%
  group_by(region, year, crop) %>%
  summarise(across(c(area, production), sum), .groups = "drop") %>%
  mutate(
    yield = production / area,
    crop = fcase(
      crop == "小麦", "Wheat",
      crop == "玉米", "Maize",
      crop == "稻谷", "Rice"
    )
  ) %>%
  drop_na() %>%
  filter(year %in% 2001:2019)

stats_country <-
  stats %>%
  group_by(year, crop) %>%
  summarise(across(c(area, production), sum), .groups = "drop") %>%
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
  mutate(region = "China")

data <- qread("data/tidied.qs", nthreads = qn) %>%
  mutate(
    crop = fifelse(str_detect(crop, "Rice"), "Rice", crop),
    AOT40 = fcase(
      crop == "Maize", AOT40_7,
      crop == "Rice", AOT40_6,
      crop == "Wheat", AOT40_7 - AOT40_1
    ),
    O3 = fcase(
      crop == "Maize",
      O3_7,
      crop == "Rice",
      O3_6,
      crop == "Wheat",
      (O3_7 * (MA - `GR&EM` + 1) - O3_1) / (MA - `GR&EM` + 0)
    )
  )

data <- data %>%
  mutate(region = fcase(
    province %in% c("辽宁省", "吉林省", "黑龙江省"), "Northeast China",
    province %in% c("上海市", "江苏省", "浙江省", "安徽省", "福建省", "江西省", "山东省", "台湾省"), "East China",
    province %in% c("北京市", "天津市", "河北省", "山西省", "内蒙古自治区"), "North China",
    province %in% c("河南省", "湖北省", "湖南省"), "Central China",
    province %in% c("广东省", "广西壮族自治区", "海南省", "香港特别行政区", "澳门特别行政区"), "South China",
    province %in% c("四川省", "贵州省", "云南省", "西藏自治区", "重庆市"), "Southwest China",
    province %in% c("陕西省", "甘肃省", "青海省", "宁夏回族自治区", "新疆维吾尔自治区"), "Northwest China"
  )) %>%
  drop_na(region)

smed_region <-
  data %>%
  group_by(year, region, crop) %>%
  summarise(across(c(contains("GOSIF"), AOT40, AOD, maxtmp, PM25, O3), weighted.mean, w = fraction, na.rm = T), .groups = "drop") %>%
  filter(year %in% 2001:2019)

smed_country <-
  data %>%
  group_by(year, crop) %>%
  summarise(across(c(contains("GOSIF"), AOT40, AOD, maxtmp, PM25, O3), weighted.mean, w = fraction, na.rm = T), .groups = "drop") %>%
  filter(year %in% 2001:2019) %>%
  mutate(region = "China")

order <- c(
  "China", "North China", "Northeast China", "East China", "South China",
  "Central China", "Southwest China", "Northwest China"
)

walk(c("Maize", "Wheat", "Rice"), function(acrop) {
  p1 <- stats_region %>%
    bind_rows(stats_country) %>%
    filter(crop == acrop) %>%
    ggplot(aes(x = year, y = yield, color = factor(region, order), alpha = factor(region, order))) +
    theme_half_open(16, font_family = "Roboto Condensed") +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank()
    ) &
    ylab("Yield")

  p2 <- smed_region %>%
    bind_rows(smed_country) %>%
    filter(crop == acrop) %>%
    ggplot(aes(x = year, y = GOSIF_sum, color = factor(region, order), alpha = factor(region, order))) +
    theme_half_open(16, font_family = "Roboto Condensed") +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank()
    ) &
    ylab("SIF")

  p3 <- smed_region %>%
    bind_rows(smed_country) %>%
    filter(crop == acrop) %>%
    ggplot(aes(x = year, y = PM25, color = factor(region, order), alpha = factor(region, order))) +
    theme_half_open(16, font_family = "Roboto Condensed") +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank()
    ) &
    ylab(TeX("PM$_{2.5}$"))

  p4 <- smed_region %>%
    bind_rows(smed_country) %>%
    filter(crop == acrop) %>%
    ggplot(aes(x = year, y = O3, color = factor(region, order), alpha = factor(region, order))) +
    theme_half_open(16, font_family = "Roboto Condensed") +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank()
    ) &
    ylab("MDA8")

  p5 <- smed_region %>%
    bind_rows(smed_country) %>%
    filter(crop == acrop) %>%
    ggplot(aes(x = year, y = AOD, color = factor(region, order), alpha = factor(region, order))) +
    theme_half_open(16, font_family = "Roboto Condensed")

  p6 <- smed_region %>%
    bind_rows(smed_country) %>%
    filter(crop == acrop) %>%
    ggplot(aes(x = year, y = AOT40, color = factor(region, order), alpha = factor(region, order))) +
    theme_half_open(16, font_family = "Roboto Condensed") &
    ylab("AOT40")

  patch <-
    wrap_plots(p1, p2, p3, p4, p5, p6, ncol = 2, guides = "collect") &
      scale_color_manual(name = NULL, values = c("black", MetBrewer::met.brewer("Juarez", n = 7))) &
      scale_alpha_manual(name = NULL, values = c(1, rep(0.2, 7))) &
      background_grid() &
      geom_line() &
      geom_point(size = 2) &
      stat_poly_eq(use_label(c("eq.label", "p.value.label")), data = . %>% filter(region == "China"), family = "Roboto Condensed")

  ggsave(str_c("figures/background_line_", acrop, ".pdf"), patch, width = 2, height = 1, scale = 7)
})
