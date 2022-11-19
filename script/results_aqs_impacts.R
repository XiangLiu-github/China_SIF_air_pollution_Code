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

# aerosol -----------------------------------------------------------------

AOD <- read_rds("data/outputs/aerosol/AOD.rds")

AOD <- cldr %>%
  lazy_dt() %>%
  inner_join(AOD) %>%
  group_by(crop, x, y, year) %>%
  summarise(AOD = mean(AOD), .groups = "drop") %>%
  inner_join(fraction) %>%
  as_tibble() %>%
  nest(fdata = -c(crop, year)) %>%
  arrange(crop, year)

aod_ctr <- read_rds("data/aod_cft.rds")

aod_ctr <- aod_ctr %>%
  mutate(aod_data = pro_map(aod_data, function(adata) {
    cldr %>%
      lazy_dt() %>%
      inner_join(adata, by = c("x", "y", "month")) %>%
      group_by(crop, year, x, y) %>%
      summarise(AOD_cft = mean(AOD_cft), .groups = "drop") %>%
      as_tibble()
  })) %>%
  unnest() %>%
  nest(aod_data = -c(crop, year))

impacts <-
  reduce(list(AOD, f1, aod_ctr), inner_join)

# adata = impacts$fdata[[71]]
# coefs = impacts$coefs[[71]]
# aod_data = impacts$aod_data[[71]]

impacts <- impacts %>%
  mutate(AOD_results = pro_pmap(list(fdata, coefs, aod_data), function(adata, coefs, aod_data) {
    coefs <- coefs %>%
      select(id, term, estimate) %>%
      pivot_wider(names_from = term, values_from = estimate) %>%
      select(contains("AOD"))

    cal_data <- adata %>%
      inner_join(aod_data, by = c("x", "y")) %>%
      mutate(AOD_cft = fifelse(AOD_cft > AOD, AOD, AOD_cft))

    rel_X <- cal_data %>%
      mutate(AOD1 = AOD_cft - AOD, AOD2 = AOD_cft^2 - AOD^2) %>%
      select(AOD1, AOD2)

    stopifnot(sum(rel_X$AOD1 > 0) == 0)

    # cal relative impact
    rel_results <- tcrossprod(as.matrix(rel_X), as.matrix(coefs)) %>%
      expm1() %>%
      rowQuantiles(probs = c(0.05, 0.5, 0.95)) %>%
      as_tibble() %>%
      bind_cols(cal_data %>% select(x, y, fraction, pm_level))

    return(rel_results)
  })) %>%
  select(-c(fdata, coefs, aod_data))

qsave(impacts, "data/impacts_aerosol.qs", nthread = qn)

impacts$AOD_results[[20]] %>%
  filter(pm_level == 50) %>%
  select(x, y, `50%`) %>%
  rast(type = "xyz") %>%
  plot()


# ozone -------------------------------------------------------------------

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

ozone <- cldr %>%
  lazy_dt() %>%
  inner_join(ozone) %>%
  group_by(crop, x, y, year) %>%
  summarise(AOT40 = sum(AOT40), .groups = "drop") %>%
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
      summarise(AOT40_cft = sum(AOT40_cft), .groups = "drop") %>%
      as_tibble()
  })) %>%
  unnest() %>%
  nest(ozone_data = -c(crop, year))

impacts <-
  reduce(list(ozone, f1, ozone_ctr), inner_join)

# adata = impacts$fdata[[1]]
# coefs = impacts$coefs[[1]]
# ozone_data = impacts$ozone_data[[1]]

impacts <- impacts %>%
  mutate(ozone_results = pro_pmap(list(fdata, coefs, ozone_data), function(adata, coefs, ozone_data) {
    coefs <- coefs %>%
      select(id, term, estimate) %>%
      pivot_wider(names_from = term, values_from = estimate) %>%
      select(contains("AOT40"))

    cal_data <- adata %>%
      inner_join(ozone_data, by = c("x", "y")) %>%
      mutate(AOT40_cft = fifelse(AOT40_cft > AOT40, AOT40, AOT40_cft))

    rel_X <- cal_data %>%
      mutate(AOT401 = AOT40_cft - AOT40) %>%
      select(AOT401)

    stopifnot(sum(rel_X$AOT401 > 0) == 0)

    # cal relative impact
    rel_results <- tcrossprod(as.matrix(rel_X), as.matrix(coefs)) %>%
      expm1() %>%
      rowQuantiles(probs = c(0.05, 0.5, 0.95)) %>%
      as_tibble() %>%
      bind_cols(cal_data %>% select(x, y, fraction, peak_level))

    return(rel_results)
  })) %>%
  select(-c(fdata, coefs, ozone_data))

qsave(impacts, "data/impacts_ozone_AOT40.qs", nthread = qn)

impacts$ozone_results[[60]] %>%
  filter(peak_level == 30) %>%
  select(x, y, `50%`) %>%
  rast(type = "xyz") %>%
  plot()


# map --------------------------------------------------------------------

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
  mutate(AOD_results = pro_map(AOD_results, function(adata) {
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
  }))

ozone <- qread("data/impacts_ozone_AOT40.qs", nthreads = qn)

ozone <-
  ozone %>%
  mutate(ozone_results = pro_map(ozone_results, function(adata) {
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
  }))

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
      legend.position = c(0.2, 0.05),
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
    xlab(NULL) +
    ylab(NULL)

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
    xlab(NULL) +
    ylab(NULL)

  return(lst(p1, p2))
}

op <- pmap(list(unique(ozone$crop_parent), c(F, F, F), letters[4:6]), ozone_plot)

ap[[2]]$p1 <- ap[[2]]$p1 +
  geom_text_npc(
    aes(npcx = 0.32, npcy = 0.02, label = "Gain"),
    family = "Roboto Condensed",
    size = 5,
    color =  WrensBookshelf::WB_brewer("BabyWrenAndTheGreatGift", direction = -1)[9]
  ) +
  geom_text_npc(
    aes(npcx = 0.04, npcy = 0.02, label = "Loss"),
    family = "Roboto Condensed",
    size = 5,
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
  width = 2,
  height = 3,
  scale = 9
)

# matrix ------------------------------------------------------------------

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
  facet_grid2(vars(region), vars(crop_parent),
    strip = strip_themed(text_y = map(
      MetBrewer::met.brewer("Juarez", n = 7) %>% append("black", 3),
      ~ element_text(color = .x)
    ))
  ) +
  geom_tile() +
  geom_contour(aes(z = value), color = "grey", bins = 20) +
  geom_point(data = . %>% filter(pm_level == 35 &
    peak_level == 60)) +
  geom_text(
    aes(label = sprintf("%0.2f", round(value, digits = 2))),
    data = . %>% filter(pm_level == 35 & peak_level == 60),
    family = "Roboto Condensed",
    nudge_y = 10,
    nudge_x = 10
  ) +
  scale_fill_gradientn(
    name = TeX("Percentage change in SIF"),
    limits = c(-10, 10),
    oob = squish,
    colours = c('#762a83', '#c51b7d', '#f7f7f7', '#5aae61', '#01665e')
  ) +
  theme_half_open(16, font_family = "Roboto Condensed") +
  background_grid() +
  scale_x_continuous(
    name = TeX("PM$_{2.5}$ level ($\\mu g$ $m^{-3}$)"),
    expand = expansion()
  ) +
  scale_y_continuous(
    name = TeX("Peak season ozone level ($\\mu$g m$^{-3}$)"),
    expand = expansion()
  ) +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(2.5, "line"),
    strip.background = element_blank()
  )

ggsave(
  "figures/hist_impact_matrix.pdf",
  p3,
  width = 3,
  height = 8,
  scale = 2
)




