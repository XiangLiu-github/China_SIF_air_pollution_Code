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
  filter(pm_level == 35) %>%
  rename_with(~ str_c("aer", .x), contains("%")) %>%
  select(-pm_level)

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
  filter(peak_level == 60) %>%
  rename_with(~ str_c("ozo", .x), contains("%")) %>%
  select(-peak_level)

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

p1 <- joined %>%
  ggplot(aes(x = year)) +
  facet_wrap(~crop, scales = "free_y") +
  geom_ribbon(aes(
    ymax = yield + `aer95%` * yield + `ozo95%` * yield,
    ymin = yield + `aer5%` * yield + `ozo5%` * yield
  ),
  alpha = 0.5, show.legend = F, color = NA, fill = PNWColors::pnw_palette("Shuksan2", n = 4)[1]
  ) +
  geom_ribbon(aes(
    ymax = yield + `aer95%` * yield,
    ymin = yield + `aer5%` * yield
  ),
  alpha = 0.2, show.legend = F, color = NA, fill = PNWColors::pnw_palette("Shuksan2", n = 4)[2]
  ) +
  geom_ribbon(aes(
    ymax = yield + `ozo95%` * yield,
    ymin = yield + `ozo5%` * yield
  ),
  alpha = 0.2, show.legend = F, color = NA, fill = PNWColors::pnw_palette("Shuksan2", n = 4)[3]
  ) +
  ggborderline::geom_borderline(aes(y = yield + `aer50%` * yield + `ozo50%` * yield), size = 1, show.legend = F, linetype = "dashed", lineend = "round", color = PNWColors::pnw_palette("Shuksan2", n = 4)[1]) +
  ggborderline::geom_borderline(aes(y = yield + `aer50%` * yield), size = 1, show.legend = F, lineend = "round", bordersize = .6, color = PNWColors::pnw_palette("Shuksan2", n = 4)[2]) +
  ggborderline::geom_borderline(aes(y = yield + `ozo50%` * yield), size = 1, show.legend = F, lineend = "round", bordersize = .6, color = PNWColors::pnw_palette("Shuksan2", n = 4)[3]) +
  geom_line(aes(y = yield), show.legend = F, color = PNWColors::pnw_palette("Shuksan2", n = 4)[4]) +
  scale_y_continuous(name = "Yield (kg/ha)") +
  scale_x_continuous(name = NULL, breaks = c(2005, 2010, 2015, 2019)) +
  theme_half_open(16, font_family = "Roboto Condensed") +
  background_grid()

p2 <- joined %>%
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
    kcal_percapita_perday_our_aer = 0.44 * area * (yield + `aer50%` * yield) * 1e3 * n * (1 - w) * E / (population * 1e4) / 365,
    kcal_percapita_perday_our_ozo = 0.44 * area * (yield + `ozo50%` * yield) * 1e3 * n * (1 - w) * E / (population * 1e4) / 365,
    kcal_percapita_perday_our_up = 0.44 * area * (yield + `aer95%` * yield + `ozo95%` * yield) * 1e3 * n * (1 - w) * E / (population * 1e4) / 365,
    kcal_percapita_perday_our_lw = 0.44 * area * (yield + `aer5%` * yield + `ozo5%` * yield) * 1e3 * n * (1 - w) * E / (population * 1e4) / 365,
  ) %>%
  group_by(year) %>%
  summarise(across(contains("kcal"), sum)) %>%
  mutate(year = as.character(year)) %>%
  bind_rows((.) %>% summarise(across(everything(), mean)) %>% mutate(year = "Avrg.")) %>%
  ggplot(aes(x = year)) +
  geom_col(aes(y = kcal_percapita_perday_our), fill = PNWColors::pnw_palette("Shuksan2", n = 4)[1], width = 0.6, color = "black") +
  geom_text(aes(y = kcal_percapita_perday_our, label = str_c("Our. ", round(kcal_percapita_perday_our))),
            data = . %>% filter(year == "Avrg."),
            family = "Roboto Condensed", nudge_x = 1.2,
            color = PNWColors::pnw_palette("Shuksan2", n = 4)[1]
  ) +
  geom_errorbar(aes(ymax = kcal_percapita_perday_our_up, ymin = kcal_percapita_perday_our_lw), width = 0.3) +
  # geom_col(aes(y = kcal_percapita_perday_our_aer), fill = PNWColors::pnw_palette("Shuksan2", n = 4)[2], width = 0.6, color = "black") +
  # geom_col(aes(y = kcal_percapita_perday_our_ozo), fill = PNWColors::pnw_palette("Shuksan2", n = 4)[3], width = 0.6, color = "black") +
  geom_col(aes(y = kcal_percapita_perday), fill = PNWColors::pnw_palette("Shuksan2", n = 4)[4], width = 0.6, color = "black") +
  geom_text(aes(y = kcal_percapita_perday, label = str_c("Obs. ", round(kcal_percapita_perday))),
            data = . %>% filter(year == "Avrg."),
            family = "Roboto Condensed", nudge_x = 1.2,
            color = PNWColors::pnw_palette("Shuksan2", n = 4)[4]
  ) +
  scale_y_continuous(name = TeX("Kcal capita$^{-1}$ day$^{-1}$")) +
  scale_x_discrete(name = NULL, expand = expansion(add = c(0.6, 2))) +
  coord_cartesian(ylim = c(850, NA)) +
  theme_half_open(16, font_family = "Roboto Condensed") +
  background_grid()

patch <- wrap_plots(p1, p2,
                    design =
                      "A
                   B
                   B"
) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 20))

ggsave("figures/hist_impact_calo_new.pdf", patch, width = 3, height = 2, scale = 3)
