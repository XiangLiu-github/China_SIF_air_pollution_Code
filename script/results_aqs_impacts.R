source("script/loadPackages.R")
source("script/loadFunctions.R")

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

aerosol <- qread("data/impacts_aerosol.qs", nthreads = qn)

aerosol <-
  aerosol %>%
  filter(year %in% 2005:2019) %>%
  mutate(AOD_results = map(AOD_results, function(adata) {
    temp <- adata %>%
      pivot_wider(names_from = pm_level, values_from = contains("%")) %>%
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
      pivot_longer(
        contains("%"),
        names_prefix = "weighted_mean.",
        names_to = c("name", "pm_level"),
        names_sep = "_"
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
      pivot_longer(
        contains("%"),
        names_prefix = "weighted_mean.",
        names_to = c("name", "pm_level"),
        names_sep = "_"
      ) %>%
      pivot_wider()

    lst(province_level, region_level)
  }, .progress = T))

ozone <- qread("data/impacts_ozone_AOT40.qs", nthreads = qn)

ozone <-
  ozone %>%
  mutate(ozone_results = map(ozone_results, function(adata) {
    temp <- adata %>%
      pivot_wider(names_from = peak_level, values_from = contains("%")) %>%
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
      pivot_longer(
        contains("%"),
        names_prefix = "weighted_mean.",
        names_to = c("name", "peak_level"),
        names_sep = "_"
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
      pivot_longer(
        contains("%"),
        names_prefix = "weighted_mean.",
        names_to = c("name", "peak_level"),
        names_sep = "_"
      ) %>%
      pivot_wider()

    lst(province_level, region_level)
  }, .progress = T))

saveRDS(aerosol, "data/impacts_aerosol_summarised.rds")
saveRDS(ozone, "data/impacts_ozone_summarised.rds")

aerosol <- read_rds("data/impacts_aerosol_summarised.rds")
ozone <- read_rds("data/impacts_ozone_summarised.rds")

aerosol_plot <- function(acrop, plot_legend, subtitle) {
  p1 <- aerosol %>%
    filter(crop_parent == acrop) %>%
    unnest_wider(AOD_results) %>%
    select(crop, year, crop_parent, province_level) %>%
    unnest() %>%
    group_by(crop_parent, 省, pm_level) %>%
    summarise(across(contains("%"), fmean), .groups = "drop") %>%
    filter(pm_level == 35) %>%
    inner_join(province_plot, c("省" = "name")) %>%
    ggplot() +
    geom_sf(
      aes(fill = `50%` * 1e2, geometry = geometry),
      size = 0.2,
      show.legend = plot_legend,
      color = NA
    ) +
    geom_sf(
      aes(color = ""),
      data = region_plot,
      fill = NA,
      size = 0.8,
      show.legend = plot_legend
    ) +
    scale_fill_gradientn(
      name = TeX("Percentage change in SIF"),
      limits = c(-10, 10),
      oob = squish,
      colours = WrensBookshelf::WB_brewer("BabyWrenAndTheGreatGift", direction = -1),
      na.value = "grey90",
      guide = guide_colorbar(
        title.position = "top",
        barwidth = 15,
        label.position = "bottom",
        order = 1
      )
    ) +
    xlab(NULL) +
    ylab(NULL) +
    theme_ipsum_rc(
      panel_spacing = grid::unit(0, "lines"),
      plot_margin = margin(0, 0, 0, 0),
      grid = F,
      axis = F,
      base_size = 15,
      axis_title_size = 15,
      strip_text_size = 20,
    ) +
    theme(
      legend.position = c(0.2, 0.1),
      legend.direction = "horizontal",
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      plot.tag = element_text(size = 50),
      legend.title.align = 0.5
    ) +
    scale_colour_manual(
      values = NA,
      guide = guide_legend("No data", order = 2)
    ) +
    labs(tag = subtitle) +
    ggnewscale::new_scale_color() +
    geom_sf(
      aes(color = factor(region, order)),
      inherit.aes = F,
      show.legend = F,
      data = region_plot %>% st_centroid() %>% filter(!is.na(region)),
      stroke = 2,
      fill = "white",
      shape = 21,
      size = 2
    ) +
    scale_color_manual(
      values = MetBrewer::met.brewer("Juarez", n = 7),
      guide = NULL
    )

  p2 <- aerosol %>%
    filter(crop_parent == acrop) %>%
    unnest_wider(AOD_results) %>%
    select(crop, year, crop_parent, region_level) %>%
    unnest() %>%
    drop_na() %>%
    group_by(crop_parent, region, pm_level) %>%
    summarise(across(contains("%"), fmean), .groups = "drop") %>%
    mutate(
      pm_level = as.integer(pm_level),
      region = factor(region, order),
      across(contains("%"), ~ .x * 1e2)
    ) %>%
    filter(pm_level <= 50) %>%
    ggplot(aes(
      x = pm_level,
      y = `50%`,
      color = region
    )) +
    facet_manual(
      vars(region),
      design = design,
      scales = "free_y",
      strip = strip_themed(text_x = map(
        MetBrewer::met.brewer("Juarez", n = 7) %>% append("black", 3),
        ~ element_text(color = .x)
      ))
    ) +
    geom_vline(xintercept = 35, linetype = "dashed") +
    geom_interval(
      aes(ymin = `5%`, ymax = `95%`),
      alpha = 0.5,
      size = 2,
      show.legend = F
    ) +
    geom_line(show.legend = F) +
    geom_point(show.legend = F, color = "black") +
    scale_color_manual(values = MetBrewer::met.brewer("Juarez", n = 7) %>% append("black", 3)) +
    theme_half_open(16, font_family = "Roboto Condensed") +
    background_grid() +
    theme(strip.background = element_blank()) +
    xlab(TeX("Annual PM$_{2.5}$ ($\\mu g$ $m^{-3}$)")) +
    ylab("Percentage change in SIF") +
    theme(axis.title = element_text(size = 25))

  return(lst(p1, p2))
}

ap <- pmap(list(unique(aerosol$crop_parent), c(F, T, F), letters[1:3]), aerosol_plot)

ozone_plot <- function(acrop, plot_legend, subtitle) {
  p1 <- ozone %>%
    filter(crop_parent == acrop) %>%
    unnest_wider(ozone_results) %>%
    select(crop, year, crop_parent, province_level) %>%
    unnest() %>%
    group_by(crop_parent, 省, peak_level) %>%
    summarise(across(contains("%"), fmean), .groups = "drop") %>%
    filter(peak_level == 60) %>%
    inner_join(province_plot, c("省" = "name")) %>%
    ggplot() +
    geom_sf(
      aes(fill = `50%` * 1e2, geometry = geometry),
      size = 0.2,
      show.legend = plot_legend,
      color = NA
    ) +
    geom_sf(
      aes(color = ""),
      data = region_plot,
      fill = NA,
      size = 0.8,
      show.legend = plot_legend
    ) +
    scale_fill_gradientn(
      name = TeX("percentage change in SIF"),
      limits = c(-10, 10),
      oob = squish,
      colours = WrensBookshelf::WB_brewer("BabyWrenAndTheGreatGift", direction = -1),
      na.value = "grey90",
      guide = guide_colorbar(
        title.position = "top",
        barwidth = 10,
        label.position = "bottom",
        order = 1
      )
    ) +
    xlab(NULL) +
    ylab(NULL) +
    theme_ipsum_rc(
      panel_spacing = grid::unit(0, "lines"),
      plot_margin = margin(0, 0, 0, 0),
      grid = F,
      axis = F,
      base_size = 15,
      axis_title_size = 15,
      strip_text_size = 20,
    ) +
    theme(
      legend.position = c(0.1, 0.1),
      legend.direction = "horizontal",
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      plot.tag = element_text(size = 50)
    ) +
    scale_colour_manual(
      values = NA,
      guide = guide_legend("No data", order = 2)
    ) +
    labs(tag = subtitle) +
    ggnewscale::new_scale_color() +
    geom_sf(
      aes(color = factor(region, order)),
      inherit.aes = F,
      show.legend = F,
      data = region_plot %>% st_centroid() %>% filter(!is.na(region)),
      stroke = 2,
      fill = "white",
      shape = 21,
      size = 2
    ) +
    scale_color_manual(
      values = MetBrewer::met.brewer("Juarez", n = 7),
      guide = NULL
    )

  p2 <- ozone %>%
    filter(crop_parent == acrop) %>%
    unnest_wider(ozone_results) %>%
    select(crop, year, crop_parent, region_level) %>%
    unnest() %>%
    drop_na() %>%
    group_by(crop_parent, region, peak_level) %>%
    summarise(across(contains("%"), fmean), .groups = "drop") %>%
    mutate(
      peak_level = as.integer(peak_level),
      region = factor(region, order),
      across(contains("%"), ~ .x * 1e2)
    ) %>%
    filter(peak_level <= 80) %>%
    ggplot(aes(
      x = peak_level,
      y = `50%`,
      color = region
    )) +
    facet_manual(
      vars(region),
      design = design,
      scales = "free_y",
      strip = strip_themed(text_x = map(
        MetBrewer::met.brewer("Juarez", n = 7) %>% append("black", 3),
        ~ element_text(color = .x)
      ))
    ) +
    geom_vline(xintercept = 60, linetype = "dashed") +
    geom_interval(
      aes(ymin = `5%`, ymax = `95%`),
      alpha = 0.5,
      size = 2,
      show.legend = F
    ) +
    geom_line(show.legend = F) +
    geom_point(show.legend = F, color = "black") +
    scale_color_manual(values = MetBrewer::met.brewer("Juarez", n = 7) %>% append("black", 3)) +
    theme_half_open(16, font_family = "Roboto Condensed") +
    background_grid() +
    theme(strip.background = element_blank()) +
    xlab(TeX("Peak season MDA8 ($\\mu$g m$^{-3}$)")) +
    ylab("Percentage change in SIF") +
    theme(axis.title = element_text(size = 25))

  return(lst(p1, p2))
}

op <- pmap(list(unique(ozone$crop_parent), c(F, F, F), letters[4:6]), ozone_plot)

ap[[2]]$p1 <- ap[[2]]$p1 +
  geom_text_npc(
    aes(npcx = 0.28, npcy = 0.04, label = "Gain"),
    family = "Roboto Condensed",
    size = 10,
    color = WrensBookshelf::WB_brewer("BabyWrenAndTheGreatGift", direction = -1)[9]
  ) +
  geom_text_npc(
    aes(npcx = 0.04, npcy = 0.04, label = "Loss"),
    family = "Roboto Condensed",
    size = 10,
    color = WrensBookshelf::WB_brewer("BabyWrenAndTheGreatGift", direction = -1)[1]
  )

patch <- map(
  c(ap, op),
  ~ .x$p1 + inset_element(
    .x$p2,
    left = 0,
    bottom = 0,
    right = 1,
    top = 1
  )
) %>%
  wrap_plots(ncol = 2, byrow = F)

ggsave(
  "figures/hist_impact_map.pdf",
  patch,
  width = 2.1,
  height = 3,
  scale = 9
)
