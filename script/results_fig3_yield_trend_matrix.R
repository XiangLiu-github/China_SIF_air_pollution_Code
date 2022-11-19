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

ozone <- qread("data/impacts_ozone_AOT40.qs", nthreads = qn)

ozone <- ozone %>%
  mutate(ozone_results = pro_map(ozone_results, function(adata) {
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

# p1 ----------------------------------------------------------------------

p1 <- joined %>%
  filter(pm_level == 35, peak_level == 60) %>% 
  mutate(obs = yield,
         counter = yield + `aer50%` * yield + `ozo50%` * yield,
         # counter_5 = yield + `aer5%` * yield + `ozo5%` * yield,
         # counter_95 = yield + `aer95%` * yield + `ozo95%` * yield,
         aer = yield + `aer50%` * yield,
         # aer_5 = yield + `aer5%` * yield,
         # aer_95 = yield + `aer95%` * yield,
         ozo = yield + `ozo50%` * yield,
         # ozo_5 = yield + `ozo5%` * yield,
         # ozo_95 = yield + `ozo95%` * yield
         ) %>% 
  pivot_longer(c(obs, counter, aer, ozo)) %>% 
  mutate(name = factor(name, c('obs', 'counter', 'aer', 'ozo'))) %>% 
  ggplot(aes(x = year, y = value, color = name)) +
  facet_wrap(~crop, scales = "free_y") +
  geom_line(size = 1) +
  scale_color_manual(values = PNWColors::pnw_palette("Shuksan2", n = 4),
                     name = NULL, labels = c('Obs.', 'A + O', 'Aerosol', 'Ozone')) +
  scale_y_continuous(name = "Yield (kg/ha)") +
  scale_x_continuous(name = NULL, breaks = c(2005, 2010, 2015, 2019)) +
  theme_half_open(16, font_family = "Roboto Condensed") +
  background_grid(major = 'x') +
  theme(strip.background = element_blank())


# p2 ----------------------------------------------------------------------

p2 = joined %>% 
  group_by(crop, pm_level, peak_level) %>% 
  summarise(across(contains('50'), mean), .groups = 'drop') %>% 
  mutate(value = (`aer50%` + `ozo50%`) * 1e2) %>% 
  ggplot(aes(x = pm_level, y = peak_level, fill = value)) +
  facet_wrap(~crop) +
  geom_tile() +
  geom_contour(aes(z = value), color = "grey", bins = 20) +
  geom_point(data = . %>% filter(pm_level == 35 &
                                   peak_level == 60)) +
  geom_star(data = . %>% group_by(crop) %>% slice_max(value, n = 1) %>% ungroup(), fill = 'yellow', color = NA, size = 2) +
  geom_text(
    aes(label = sprintf("%0.2f", round(value, digits = 2))),
    data = . %>% filter(pm_level == 35 & peak_level == 60),
    size = 5,
    family = "Roboto Condensed",
    nudge_y = 10,
    nudge_x = 10
  ) +
  geom_path(aes(x = PM25, y = O3, alpha = year), color = 'red',
            data = hist_point %>% filter(region == 'China') %>% rename(crop = crop_parent), show.legend = F, inherit.aes = F) +
  geom_text_repel(aes(x = PM25, y = O3, label = year), color = 'red',
            data = hist_point %>% filter(region == 'China', year %in% c(2005, 2019)) %>% rename(crop = crop_parent),
            inherit.aes = F, family = "Roboto Condensed", box.padding = 0.5, min.segment.length = Inf,
            direction = 'both', vjust = 'top', hjust = 'left') +
  scale_fill_gradientn(
    name = TeX("Percentage change in yield"),
    limits = c(-10, 10),
    oob = squish,
    colours = c('#762a83', '#c51b7d', '#f7f7f7', '#5aae61', '#01665e')
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
    legend.justification = 'center',
    legend.title = element_text(angle = -90, family = "Roboto Condensed"),
    strip.text = element_blank(), 
    legend.key.height = unit(2.5, 'lines'),
    plot.tag.position = 'topleft'
  ) +
  guides(fill = guide_colorbar(title.position = "right"))

patch = p1 / p2 +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 20))

ggsave("figures/hist_impact_fig3.pdf", patch, width = 3.5, height = 2.5, scale = 3)


