source("script/loadPackages.R")
source("script/loadFunctions.R")

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

pm <- read_rds("data/outputs/aerosol/PM25.rds") %>%
  fgroup_by(x, y, year) %>%
  fsummarise(PM25 = fmean(PM25)) %>%
  pivot_wider(names_from = year, values_from = PM25) %>%
  rast(type = "xyz", crs = "epsg:4326")

pm <-
  tibble(
    crop_parent = c("Maize", "Rice", "Wheat"),
    data = map(crop_parent, function(acrop) {
      amask <- rast(str_c("data/outputs/masks/mask_", acrop, ".tif"))

      exact_extract(pm[[1:20]], region, "weighted_mean",
        stack_apply = T,
        append_cols = "region",
        weights = amask, default_weight = 0, progress = F
      ) %>%
        pivot_longer(-region,
          names_to = c("year", NA), names_prefix = "weighted_mean.", names_sep = "\\.", names_transform = list(year = as.integer),
          values_to = "PM25"
        )
    })
  ) %>%
  unnest()

ozo <- read_rds("data/outputs/ozone/tidied.rds") # unit is ug/m3

ozo <- ozo %>%
  fsubset(month %in% 4:9) %>%
  fgroup_by(x, y, year) %>%
  fsummarise(O3 = fmean(O3)) %>%
  pivot_wider(names_from = year, values_from = O3) %>%
  rast(type = "xyz", crs = "epsg:4326")

ozo <-
  tibble(
    crop_parent = c("Maize", "Rice", "Wheat"),
    data = map(crop_parent, function(acrop) {
      amask <- rast(str_c("data/outputs/masks/mask_", acrop, ".tif"))

      exact_extract(ozo[[1:15]], region, "weighted_mean",
        stack_apply = T,
        append_cols = "region",
        weights = amask[[6:20]], default_weight = 0, progress = F
      ) %>%
        pivot_longer(-region,
          names_to = c("year", NA), names_prefix = "weighted_mean.", names_sep = "\\.", names_transform = list(year = as.integer),
          values_to = "O3"
        )
    })
  ) %>%
  unnest()

hist_point <- inner_join(ozo, pm) %>%
  mutate(region = factor(region, order))

# first run results_aqs_impacts.R
aerosol <- read_rds("data/impacts_aerosol_hete_full_summarised.rds")
ozone <- read_rds("data/impacts_ozone_hete_full_summarised.rds")

design <-
  "
ABCMNO
DEFPQR
GHISTU
JKLVWX
"

p3 <- inner_join(
  aerosol %>%
    unnest_wider(AOD_results) %>%
    select(crop, year, crop_parent, region_level) %>%
    unnest() %>%
    drop_na() %>%
    group_by(crop_parent, region, pm_level) %>%
    summarise(across(contains("%"), fmean), .groups = "drop") %>%
    mutate(
      region = factor(region, order),
      across(contains("%"), ~ .x * 1e2)
    ) %>%
    pivot_longer(contains("%"), values_to = "aerosol"),
  ozone %>%
    unnest_wider(ozone_results) %>%
    select(crop, year, crop_parent, region_level) %>%
    unnest() %>%
    drop_na() %>%
    group_by(crop_parent, region, peak_level) %>%
    summarise(across(contains("%"), fmean), .groups = "drop") %>%
    mutate(
      region = factor(region, order),
      across(contains("%"), ~ .x * 1e2)
    ) %>%
    pivot_longer(contains("%"), values_to = "ozone")
) %>%
  filter(name == "50%") %>%
  mutate(
    value = aerosol + ozone,
    across(contains("level"), as.integer)
  ) %>%
  ggplot(aes(x = pm_level, y = peak_level, fill = value)) +
  facet_manual(vars(region, crop_parent),
    design = design,
    strip = strip_themed(text_x = map(
      rep(MetBrewer::met.brewer("Juarez", n = 7) %>% append("black", 3), each = 3),
      ~ element_text(color = .x)
    ))
  ) +
  geom_tile() +
  geom_contour(aes(z = value), color = "grey", bins = 20) +
  geom_point(data = . %>% filter(pm_level == 35 &
    peak_level == 60)) +
  geom_star(data = . %>% group_by(crop_parent, region) %>% slice_max(value, n = 1) %>% ungroup(), fill = "yellow", color = NA, size = 2) +
  geom_text(
    aes(label = sprintf("%0.2f", round(value, digits = 2))),
    data = . %>% filter(pm_level == 35 & peak_level == 60),
    family = "Roboto Condensed",
    nudge_y = 10,
    nudge_x = 10
  ) +
  geom_text_repel(
    aes(label = sprintf("%0.2f", round(value, digits = 2))),
    data = . %>% slice_max(value, n = 1, by = c(crop_parent, region)), color = "yellow",
    family = "Roboto Condensed"
  ) +
  geom_path(aes(x = PM25, y = O3, alpha = year, group = 1),
    color = "red",
    data = hist_point, show.legend = F,
    inherit.aes = F
  ) +
  geom_text_repel(aes(x = PM25, y = O3, label = year),
    color = "red",
    data = hist_point %>% filter(year %in% c(2005, 2019)),
    inherit.aes = F, family = "Roboto Condensed",
    box.padding = 0.5, min.segment.length = Inf, size = 4,
    direction = "both", vjust = "top", hjust = "left"
  ) +
  scale_fill_gradientn(
    name = TeX("Percentage change in SIF"),
    limits = c(-10, 10),
    oob = squish,
    colours = c("#762a83", "#c51b7d", "#f7f7f7", "#5aae61", "#01665e")
  ) +
  theme_half_open(16, font_family = "Roboto Condensed") +
  background_grid() +
  scale_x_continuous(
    name = TeX("Annual PM$_{2.5}$ ($\\mu g$ $m^{-3}$)"),
    expand = expansion()
  ) +
  scale_y_continuous(
    name = TeX("Peak season MDA8 ($\\mu$g m$^{-3}$)"),
    expand = expansion()
  ) +
  theme(
    legend.position = "bottom", legend.justification = "center",
    legend.key.width = unit(2.5, "line"),
    strip.background = element_blank()
  )

ggsave(
  "figures/hist_impact_matrix_new.pdf",
  p3,
  width = 6,
  height = 6,
  scale = 2
)
