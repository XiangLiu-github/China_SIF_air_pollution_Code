source("script/loadPackages.R")
source("script/loadFunctions.R")

province_plot <- st_read("https://geo.datav.aliyun.com/areas_v3/bound/100000_full.json") %>%
  mutate(region = fcase(
    name %in% c("辽宁省", "吉林省", "黑龙江省"), "Northeast China",
    name %in% c("上海市", "江苏省", "浙江省", "安徽省", "福建省", "江西省", "山东省", "台湾省"), "East China",
    name %in% c("北京市", "天津市", "河北省", "山西省", "内蒙古自治区"), "North China",
    name %in% c("河南省", "湖北省", "湖南省"), "Central China",
    name %in% c("广东省", "广西壮族自治区", "海南省", "香港特别行政区", "澳门特别行政区"), "South China",
    name %in% c("四川省", "贵州省", "云南省", "西藏自治区", "重庆市"), "Southwest China",
    name %in% c("陕西省", "甘肃省", "青海省", "宁夏回族自治区", "新疆维吾尔自治区"), "Northwest China"
  ))

aod <- read_rds("data/outputs/aerosol/AOD.rds") %>%
  fgroup_by(year, x, y) %>%
  fsummarise(AOD = fmean(AOD)) %>%
  filter(year %in% 2005:2019) %>%
  fgroup_by(x, y) %>%
  fsummarise(`Observed level` = fmean(AOD))

aod_cft <- read_rds("data/aod_cft.rds") %>%
  filter(pm_level == 35) %>%
  pull(aod_data) %>%
  .[[1]] %>%
  summarise(`Counterfactual level` = fmean(AOD_cft), .by = c(x, y))

ozone <- read_rds("data/outputs/ozone/tidied.rds") %>%
  fgroup_by(year, x, y) %>%
  fsummarise(AOT40 = fsum(AOT40)) %>%
  filter(year %in% 2005:2019) %>%
  fgroup_by(x, y) %>%
  fsummarise(`Observed level` = fmean(AOT40))

ozone_cft <- read_rds("data/ozone_cft.rds") %>%
  filter(peak_level == 60) %>%
  pull(ozone_data) %>%
  .[[1]] %>%
  summarise(`Counterfactual level` = fsum(AOT40_cft), .by = c(x, y))

p1 <- inner_join(aod, aod_cft) %>%
  pivot_longer(-c(x, y)) %>%
  ggplot(aes(x = x, y = y, fill = value)) +
  facet_wrap(~name) +
  geom_tile() +
  geom_sf(data = province_plot %>% select(-name), inherit.aes = F, fill = NA, size = 0.8) +
  scale_fill_gradientn(
    colors = wesanderson::wes_palette("Zissou1", 10, type = "continuous"),
    name = "AOD",
    guide = guide_colorbar(
      title.position = "top", barwidth = 20,
      label.position = "bottom"
    )
  ) +
  theme_ipsum_rc(
    panel_spacing = grid::unit(0, "lines"),
    plot_margin = margin(0, 0, 0, 0),
    grid = F, axis = F,
    base_size = 15,
    axis_title_size = 15,
    strip_text_size = 20,
  ) +
  theme(
    legend.position = c(0.5, 0.2),
    legend.direction = "horizontal",
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )


p2 <- inner_join(ozone, ozone_cft) %>%
  pivot_longer(-c(x, y)) %>%
  ggplot(aes(x = x, y = y, fill = value)) +
  facet_wrap(~name) +
  geom_tile() +
  geom_sf(data = province_plot %>% select(-name), inherit.aes = F, fill = NA, size = 0.8) +
  scale_fill_gradientn(
    colors = MetBrewer::met.brewer("Hiroshige", n = 10, direction = -1),
    name = "AOT40 (ppb h)",
    limits = c(0, 6e4),
    oob = squish,
    guide = guide_colorbar(
      title.position = "top", barwidth = 20,
      label.position = "bottom"
    )
  ) +
  theme_ipsum_rc(
    panel_spacing = grid::unit(0, "lines"),
    plot_margin = margin(0, 0, 0, 0),
    grid = F, axis = F,
    base_size = 15,
    axis_title_size = 15,
    strip_text_size = 20,
  ) +
  theme(
    legend.position = c(0.5, 0.2),
    legend.direction = "horizontal",
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )

patch <- p1 / p2 +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 30))

ggsave("figures/counterfactual_map.pdf", patch, width = 1, height = 1, scale = 15)
