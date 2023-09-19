source("script/loadPackages.R")
source("script/loadFunctions.R")

# unit is the same
# luancheng c(114.683333, 37.883333)
# yucheng c(116.5702, 36.8290)

points <-
  tribble(
    ~longitude, ~latitude,
    114.4127778, 37.5313889,
    116.5702, 36.829
  ) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = "+proj=longlat +datum=WGS84 +no_defs")

# points = matrix(c(c(114.683333, 116.5702),
#                   c(37.883333, 36.8290)), nrow = 2) %>%
#   st_multipoint() %>%
#   as("Spatial") %>%
#   as("SpatialPoints") %>%
#   st_as_sf() %>%
#   `st_crs<-`("+proj=longlat +datum=WGS84 +no_defs") %>%
#   as("Spatial")

luancheng <- read_xlsx("data/inputs/EC/栾城日通量数据2007-2018.xlsx",
  sheet = "碳通量",
  range = "A1:C5000"
) %>%
  filter(is.na(...3)) %>%
  select(-...3) %>%
  rename(NEE = `NEE(gC/m2/d)`) %>%
  mutate(Date = ym(format(Date, "%Y-%m"))) %>%
  group_by(Date) %>%
  summarise(NEE = mean(NEE), .groups = "drop")

yucheng <- list.files("data/inputs/EC/禹城/通量数据/日尺度/", full.names = T) %>%
  map_dfr(read_xlsx) %>%
  filter(年 != "/") %>%
  type_convert() %>%
  mutate(Date = ymd(str_c(年, 月, 日, sep = "_"))) %>%
  select(Date, NEE, GEE) %>%
  mutate(Date = ym(format(Date, "%Y-%m"))) %>%
  group_by(Date) %>%
  summarise(NEE = mean(NEE), GEE = mean(GEE), .groups = "drop")

GOSIF <- read_rds("data/outputs/SIF/GOSIF.rds") %>%
  pivot_wider(names_from = c(year, month), values_from = GOSIF) %>%
  rast(type = "xyz", crs = "+proj=longlat +datum=WGS84 +no_defs") %>%
  terra::extract(vect(points)) %>%
  mutate(ID = as.character(ID)) %>%
  pivot_longer(-ID, names_to = "Date", names_transform = list(Date = ym), values_to = "GOSIF")

CSIF <- read_rds("data/outputs/SIF/CSIF.rds") %>%
  pivot_wider(names_from = c(year, month), values_from = CSIF) %>%
  rast(type = "xyz", crs = "+proj=longlat +datum=WGS84 +no_defs") %>%
  terra::extract(vect(points)) %>%
  mutate(ID = as.character(ID)) %>%
  pivot_longer(-ID, names_to = "Date", names_transform = list(Date = ym), values_to = "CSIF")

RTSIF <- read_rds("data/outputs/SIF/RTSIF.rds") %>%
  pivot_wider(names_from = c(year, month), values_from = RTSIF) %>%
  rast(type = "xyz", crs = "+proj=longlat +datum=WGS84 +no_defs") %>%
  terra::extract(vect(points)) %>%
  mutate(ID = as.character(ID)) %>%
  pivot_longer(-ID, names_to = "Date", names_transform = list(Date = ym), values_to = "RTSIF")

EC <- bind_rows(luancheng, yucheng, .id = "ID")

cldr <- read_rds("data/outputs/calendar/tidied.rds")

MA <- cldr %>%
  filter(crop %in% c("Maize", "Wheat")) %>%
  distinct(crop, x, y, year, MA) %>%
  pivot_wider(names_from = c(crop, year), values_from = MA) %>%
  rast(type = "xyz", crs = "+proj=longlat +datum=WGS84 +no_defs") %>%
  terra::extract(vect(points)) %>%
  mutate(ID = as.character(ID)) %>%
  pivot_longer(-ID,
    names_to = c("crop", "year"), names_sep = "_",
    values_to = "MA"
  )

GE <- cldr %>%
  filter(crop %in% c("Maize", "Wheat")) %>%
  distinct(crop, x, y, year, `GR&EM`) %>%
  pivot_wider(names_from = c(crop, year), values_from = `GR&EM`) %>%
  rast(type = "xyz", crs = "+proj=longlat +datum=WGS84 +no_defs") %>%
  terra::extract(vect(points)) %>%
  mutate(ID = as.character(ID)) %>%
  pivot_longer(-ID,
    names_to = c("crop", "year"), names_sep = "_",
    values_to = "GE"
  )

growing_season <-
  inner_join(MA, GE) %>%
  drop_na() %>%
  mutate(across(c(MA, GE), ~ ym(str_c(year, "_", .x)))) %>%
  filter(ID == "2" & year %in% 2003:2010)

joined <- EC %>%
  left_join(GOSIF) %>%
  left_join(CSIF) %>%
  left_join(RTSIF) %>%
  mutate(
    GOSIF = GOSIF / 1e4,
    GEE = -GEE
  )

scale <- 30

order <- c("GPP", "GOSIF", "RTSIF", "CSIF")

p1 <- joined %>%
  filter(ID == "2") %>%
  pivot_longer(c(GOSIF, CSIF, RTSIF), names_to = "SIF") %>%
  mutate(GEE = fifelse(GEE <= 0, 0, GEE)) %>%
  ggplot(aes(x = Date)) +
  geom_rect(aes(xmin = GE, xmax = MA, ymin = -Inf, ymax = Inf, fill = crop),
    data = growing_season, inherit.aes = F, alpha = 0.2
  ) +
  geom_line(aes(y = value, color = SIF)) +
  geom_point(aes(y = value, color = SIF)) +
  geom_line(aes(y = GEE / scale, color = "GPP")) +
  geom_point(aes(y = GEE / scale, color = "GPP")) +
  scale_y_continuous(
    name = TeX("SIF (mW $m^{-2}$ $nm^{-1}$ $sr^{-1}$)"),
    limits = c(-0.01, 0.91),
    sec.axis = sec_axis(~ . * scale, name = TeX("GPP (gC $m^{-2}$ $month^{-1}$)"))
  ) +
  scale_x_date(name = NULL, date_breaks = "1 years", labels = label_date_short("%Y"), expand = expansion(mult = 0.02, add = 1)) +
  scale_color_manual(values = c("black", "#C9CACB", "#F6CC7F", "#C57F7F", "green"), breaks = order, name = NULL) +
  theme_half_open(18, font_family = "Roboto Condensed") +
  background_grid() +
  theme(legend.title = element_blank())

formula <- y ~ x - 1

p2 <- inner_join(MA, GE) %>%
  drop_na() %>%
  mutate(month = map2(GE, MA, seq)) %>%
  unnest() %>%
  mutate(Date = ym(str_c(year, "_", month))) %>%
  select(ID, crop, Date) %>%
  left_join(EC) %>%
  left_join(GOSIF) %>%
  left_join(CSIF) %>%
  left_join(RTSIF) %>%
  mutate(
    GOSIF = GOSIF / 1e4,
    GEE = -GEE
  ) %>%
  filter(ID == "2") %>%
  drop_na() %>%
  pivot_longer(c(GOSIF, CSIF, RTSIF), names_to = "SIF") %>%
  mutate(SIF = factor(SIF, c("GOSIF", "RTSIF", "CSIF"))) %>%
  ggplot(aes(x = value, y = GEE, color = crop, fill = crop)) +
  facet_wrap(~SIF, scales = "free_x") +
  geom_point(show.legend = F) +
  stat_poly_line(formula = formula, show.legend = F) +
  stat_poly_eq(use_label(c("eq", "R2", "p.value")), formula = formula, show.legend = F, family = "Roboto Condensed") +
  theme_half_open(18, font_family = "Roboto Condensed") +
  background_grid() +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom", legend.justification = "center"
  ) +
  scale_y_continuous(limits = c(-0.01, 0.91) * scale) +
  xlab(TeX("SIF (mW $m^{-2}$ $nm^{-1}$ $sr^{-1}$)")) +
  ylab(TeX("GPP (gC $m^{-2}$ $month^{-1}$)"))

patch <- p1 / p2 + plot_layout(guides = "collect") + plot_annotation(tag_levels = "a") &
  theme(legend.position = "bottom", legend.justification = "center", plot.tag = element_text(size = 30))
ggsave("figures/EC_validation.pdf", patch, width = 2, height = 1.2, scale = 7)
