source("script/loadPackages.R")
source("script/loadFunctions.R")


# stats -------------------------------------------------------------------

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

pop <-
  read_xlsx("data/inputs/province_data/population.xlsx", range = "A4:AF26") %>%
  rename(year = 时间) %>%
  pivot_longer(-year, names_to = "province", values_to = "population") %>%
  mutate(year = parse_number(year))

stats <-
  inner_join(Area, Production) %>%
  inner_join(pop) %>%
  mutate(
    production = production * 1e4,
    crop = fcase(
      crop == "小麦", "Wheat",
      crop == "玉米", "Maize",
      crop == "稻谷", "Rice"
    )
  ) %>%
  drop_na() %>%
  group_by(crop, year) %>%
  summarise(
    yield = sum(production) / sum(area),
    area = sum(area),
    population = sum(population),
    .groups = "drop"
  )


# histrical pm and ozo level ----------------------------------------------

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

pm = 
  tibble(crop_parent = c('Maize', 'Rice', 'Wheat'),
         data = pro_map(crop_parent, function(acrop){
           
           amask = rast(str_c('data/outputs/masks/mask_', acrop, '.tif'))
           
           exact_extract(pm[[1:20]], region, 'weighted_mean', stack_apply = T, 
                         append_cols = 'region',
                         weights = amask, default_weight = 0, progress = F) %>% 
             pivot_longer(-region,
                          names_to = c('year', NA), names_prefix = 'weighted_mean.', names_sep = '\\.', names_transform = list(year = as.integer),
                          values_to = 'PM25')
           
         })) %>% 
  unnest()

ozo <- read_rds("data/outputs/ozone/tidied.rds") # unit is ug/m3

ozo <- ozo %>%
  fsubset(month %in% 4:9) %>% 
  fgroup_by(x, y, year) %>% 
  fsummarise(O3 = fmean(O3)) %>% 
  pivot_wider(names_from = year, values_from = O3) %>% 
  rast(type = "xyz", crs = "epsg:4326")

ozo = 
  tibble(crop_parent = c('Maize', 'Rice', 'Wheat'),
         data = pro_map(crop_parent, function(acrop){
           
           amask = rast(str_c('data/outputs/masks/mask_', acrop, '.tif'))
           
           exact_extract(ozo[[1:15]], region, 'weighted_mean', stack_apply = T, 
                         append_cols = 'region',
                         weights = amask[[6:20]], default_weight = 0, progress = F) %>% 
             pivot_longer(-region,
                          names_to = c('year', NA), names_prefix = 'weighted_mean.', names_sep = '\\.', names_transform = list(year = as.integer),
                          values_to = 'O3')
           
         })) %>% 
  unnest()

hist_point = inner_join(ozo, pm) %>% 
  mutate(region = factor(region, order))


# counterfacutal ----------------------------------------------------------

aerosol <- qread("data/impacts_aerosol.qs", nthreads = qn)

aerosol <- aerosol %>%
  mutate(AOD_results = pro_map(AOD_results, function(adata) {
    adata %>%
      group_by(pm_level) %>%
      summarise(across(contains("%"), ~ weighted.mean(.x, fraction)))
  })) %>%
  unnest() %>%
  group_by(crop_parent, year, pm_level) %>%
  summarise(across(contains("%"), mean), .groups = "drop") %>%
  rename_with(~ str_c("aer", .x), contains("%")) 

ozone <- qread("data/impacts_ozone_Feng2022.qs", nthreads = qn)

ozone <- ozone %>%
  mutate(
    crop_parent = fifelse(str_detect(crop, "Rice"), "Rice", crop),
    ozone_results = pro_map(ozone_results, function(adata) {
    adata %>%
      group_by(peak_level) %>%
      summarise(across(contains("%"), ~ weighted.mean(.x, fraction)))
  })) %>%
  unnest() %>%
  group_by(crop_parent, year, peak_level) %>%
  summarise(across(contains("%"), mean), .groups = "drop") %>%
  rename_with(~ str_c("ozo", .x), contains("%"))

joined <-
  inner_join(
    stats,
    aerosol,
    by = c("crop" = "crop_parent", "year")
  ) %>%
  inner_join(
    ozone,
    by = c("crop" = "crop_parent", "year")
  ) 

joined = joined %>%
  mutate(
    n = fcase(
      crop == "Wheat", 0.78,
      crop == "Rice", 1,
      crop == "Maize", 0.79
    ),
    w = fcase(
      crop == "Wheat", 0.2,
      crop == "Rice", 0.1,
      crop == "Maize", 0.7
    ),
    E = fcase(
      crop == "Wheat", 3391.67,
      crop == "Rice", 3882.05,
      crop == "Maize", 3622.95
    ),
    kcal_percapita_perday = 0.44 * area * yield * 1e3 * n * (1 - w) * E / (population * 1e4) / 365,
    kcal_percapita_perday_our = 0.44 * area * (yield + `aer50%` * yield + `ozo50%` * yield) * 1e3 * n * (1 - w) * E / (population * 1e4) / 365,
    kcal_percapita_perday_our_up = 0.44 * area * (yield + `aer95%` * yield + `ozo95%` * yield) * 1e3 * n * (1 - w) * E / (population * 1e4) / 365,
    kcal_percapita_perday_our_lw = 0.44 * area * (yield + `aer5%` * yield + `ozo5%` * yield) * 1e3 * n * (1 - w) * E / (population * 1e4) / 365,
  ) 

p1 = joined %>% 
  filter(peak_level == 60, pm_level == 35) %>% 
  group_by(year) %>%
  summarise(across(contains("kcal"), sum)) %>%
  ggplot(aes(x = year)) +
  geom_col(aes(y = kcal_percapita_perday_our), fill = '#C15555',, width = 0.5, color = "black") +
  geom_errorbar(aes(ymax = kcal_percapita_perday_our_up, ymin = kcal_percapita_perday_our_lw), width = 0.25) +
  geom_col(aes(y = kcal_percapita_perday), fill = '#89AEE5', width = 0.5, color = "black") +
  scale_y_continuous(name = TeX("kCal capita$^{-1}$ day$^{-1}$")) +
  scale_x_continuous(name = NULL, breaks = c(2005, 2010, 2015, 2019)) +
  coord_cartesian(ylim = c(850, NA)) +
  theme_half_open(16, font_family = "Roboto Condensed") +
  background_grid(major = 'y')

p2 = joined %>% 
  filter(peak_level == 60, pm_level == 35) %>% 
  group_by(crop, year) %>%
  summarise(across(contains("kcal"), sum)) %>% 
  group_by(crop) %>% 
  summarise(across(everything(), mean)) %>% 
  select(crop, kcal_percapita_perday, kcal_percapita_perday_our) %>% 
  pivot_longer(-crop) %>% 
  group_by(name) %>% 
  mutate(per = value / sum(value),
         ymax = cumsum(per),
         ymin = c(0, head(per, n = -1))) %>% 
  ungroup() %>% 
  mutate(name = fifelse(str_detect(name, 'our'), 'Counterfutal scenario', 'Observed scenario')) %>% 
  ggplot(aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3, fill = crop)) +
  facet_wrap2(~name, strip = strip_themed(text_x = map(c('#C15555', '#89AEE5'), ~ element_text(colour = .x)))) +
  geom_rect() +
  geom_text(x = 3.5, aes(y=(ymax + ymin) / 2, label = sprintf("%0.0f", round(value, digits = 2))),
            size = 5, family = "Roboto Condensed") +
  geom_text_npc(aes(npcx = 0, npcy = 0, label = sprintf("%0.0f", round(value, digits = 2))),
                data = . %>% group_by(name) %>% summarise(value = sum(value)), inherit.aes = F,
                family = "Roboto Condensed", size = 5) +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  scale_fill_manual(values = NatParksPalettes::natparks.pals("Yellowstone", n = 4)[2:4]) +
  theme_void(16, base_family = "Roboto Condensed") +
  theme(legend.title = element_blank(),
        strip.text.x = element_text(size = 16))

p3 = joined %>% 
  group_by(pm_level, peak_level, year) %>%
  summarise(across(contains("kcal"), sum), .groups = 'drop') %>% 
  group_by(pm_level, peak_level) %>%
  summarise(across(contains("kcal"), mean), .groups = 'drop') %>% 
  mutate(value = (kcal_percapita_perday_our - kcal_percapita_perday) / kcal_percapita_perday * 1e2) %>% 
  
  ggplot(aes(x = pm_level, y = peak_level, fill = value)) +
  geom_tile() +
  geom_contour(aes(z = value), color = "grey", bins = 20) +
  geom_point(data = . %>% filter(pm_level == 35, peak_level == 60)) +
  geom_star(data = . %>% slice_max(value, n = 1), fill = 'yellow', color = NA, size = 3) +
  geom_text(
    aes(label = sprintf("%0.2f", round(value, digits = 2))),
    data = . %>% filter(pm_level == 35 & peak_level == 60),
    family = "Roboto Condensed", size = 5,
    nudge_y = 5,
    nudge_x = 5
  ) +
  geom_path(aes(x = PM25, y = O3, alpha = year, group = 1), color = 'red',
            data = hist_point %>% filter(region == 'China') %>% group_by(year) %>%  summarise(across(c(O3, PM25), mean)), show.legend = F,
            inherit.aes = F) +
  geom_text_repel(aes(x = PM25, y = O3, label = year), color = 'red',
                  data = hist_point %>% filter(region == 'China', year %in% c(2005, 2019)) %>% group_by(year) %>%  summarise(across(c(O3, PM25), mean)),
                  inherit.aes = F, family = "Roboto Condensed",
                  box.padding = 0.5, min.segment.length = Inf, size = 5,
                  direction = 'both', vjust = 'top', hjust = 'left') +
  scale_fill_gradientn(
    name = TeX("Percentage change in kCal capita$^{-1}$ day$^{-1}$"),
    limits = c(0, 20),
    oob = squish,
    colours = PNWColors::pnw_palette('Lake') %>% rev
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
  theme(legend.key.width = unit(3, units = 'lines')) +
  guides(fill = guide_colorbar(title.position = "bottom"))

design = 
  '
AAA
BBC'

patch = p1 + p2 + p3 +
  plot_annotation(tag_levels = "a") +
  plot_layout(design = design, guides = 'collect') &
  theme(legend.position = 'bottom',
        plot.tag = element_text(size = 30))

ggsave("figures/hist_impact_fig4_Feng2022.pdf", patch, width = 3, height = 2.5, scale = 4)














