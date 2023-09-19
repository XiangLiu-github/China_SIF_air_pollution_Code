source("script/loadPackages.R")
source("script/loadFunctions.R")

China <- read_rds("D:\\Data/Site_air_quality/China/combinded.rds") %>%
  select(-type, -year, -DOY, -month) %>%
  as.data.table()

China <- melt(China,
  id.vars = c("hour", "date"),
  variable.name = c("site"),
  value.name = c("ppb")
)[, ppb := ppb / 2] %>%
  drop_na(ppb)

China[, `:=`(year = year(date), month = month(date), day = day(date))]

China <- China %>%
  mutate(ppb = fifelse(year >= 2018 & month >= 9, ppb * 2 / 48 * 24.46, ppb * 2 / 48 * 22.41))

China <- rollingMean(China, pollutant = "ppb")

China <- China %>%
  ungroup() %>%
  select(-date)

station <- read_excel("../Vegetation_ozone/data/inputs/station_raw/站点列表-2021.01.01起.xlsx") %>%
  rename(name = 监测点编码, lon = 经度, lat = 纬度) %>%
  select(name, lon, lat) %>%
  drop_na() %>%
  rename(site = name, Latitude = lat, Longitude = lon)

China <- inner_join(China, station) %>% select(-site)

China <- China %>% as.data.table()

Veg <- China[hour %between% c(8, 19)]
Veg[, `:=`(
  AOT40 = fifelse(ppb >= 40, ppb - 40, 0),
  W126 = ppb / (1 + 4403 * exp(-126 * ppb / 1000))
)]
gc()

Veg <- Veg %>%
  # daily
  fgroup_by(Latitude, Longitude, year, month, day) %>%
  fsummarise(AOT40 = sum(AOT40), W126 = sum(W126), data_capture = length(ppb) / 12) %>%
  filter(data_capture >= 0.75) %>%
  mutate(across(c(AOT40, W126), ~ . / data_capture)) %>%
  select(-data_capture)

MDA8 <- China[, .(o3_max = max(rolling8ppb), count = .N / 24 * 100), by = list(year, month, day, Latitude, Longitude)][count > 75][, count := NULL]

set.seed(2023)
joined <- inner_join(Veg, MDA8) %>%
  drop_na() %>%
  mutate(O3 = o3_max * 2) %>%
  slice_sample(prop = 0.01) %>%
  filter(O3 <= 400 & AOT40 <= 1e3)

require(parallel)
nc <- 9 ## cluster size, set for example portability
if (detectCores() > 1) { ## no point otherwise
  cl <- makeCluster(nc)
  ## could also use makeForkCluster, but read warnings first!
} else {
  cl <- NULL
}

folds <- vfold_cv(joined)

model_aot40 <- bam(AOT40 ~ s(O3), data = joined, family = tw(), chunk.size = 5000, cluster = cl)

# cv_w126 <- pro_map_dfr(folds$splits, function(asplit) {
#   model <- bam(W126 ~ s(O3),
#     data = asplit %>% analysis(),
#     family = tw(), chunk.size = 5000, cluster = cl
#   )
#
#   tibble(
#     predicted = exp(predict(model, asplit %>% assessment(), cluster = cl)),
#     raw = asplit %>% assessment() %>% pull(W126)
#   )
# })

cv_aot40 <- map_dfr(folds$splits, function(asplit) {
  model <- bam(AOT40 ~ s(O3),
    data = asplit %>% analysis(),
    family = tw(), chunk.size = 5000, cluster = cl
  )

  tibble(
    predicted = exp(predict(model, asplit %>% assessment(), cluster = cl)),
    raw = asplit %>% assessment() %>% pull(AOT40)
  )
}, .progress = T)


p3 <- joined %>%
  ggplot(aes(x = O3, y = AOT40)) +
  geom_pointdensity(adjust = 4, method = "default") +
  scale_color_viridis(
    trans = "log10",
    name = "n_neighbors",
    limits = c(1, 4e2), oob = scales::squish,
    option = "magma"
  ) +
  xlab(TeX("MDA8 ($\\mu g $ $ m ^{-3}$)")) +
  ylab(TeX("AOT40 (ppb h)")) +
  # scale_y_continuous(limits = c(0, 23e3)) +
  # scale_x_continuous(limits = c(0, 230)) +
  geom_xsidehistogram(bins = 60) +
  geom_ysidehistogram(bins = 60) +
  ggside(x.pos = "top", y.pos = "right", collapse = "all") +
  scale_xsidey_continuous(labels = NULL) +
  scale_ysidex_continuous(labels = NULL) +
  geom_line(aes(x = O3, y = AOT40_pred),
    inherit.aes = F,
    data = tibble(O3 = seq(0, 300, by = 0.1)) %>%
      mutate(AOT40_pred = exp(predict(model_aot40, .))),
    size = 1
  ) +
  theme_half_open(18, font_family = "Roboto Condensed") +
  background_grid() +
  theme(legend.position = c(0.2, 0.6))

p4 <- cv_aot40 %>%
  ggplot(aes(x = raw, y = predicted)) +
  geom_pointdensity(adjust = 4, method = "default") +
  scale_color_viridis(
    trans = "log10",
    name = "n_neighbors",
    limits = c(1, 1e2), oob = scales::squish,
    option = "magma"
  ) +
  geom_abline(slope = 1, intercept = 0, size = 1, color = "blue") +
  geom_smooth(method = "lm", se = F, color = "red", size = 1) +
  xlab(TeX("Observation AOT40 (ppb h)")) +
  ylab(TeX("Prediction AOT40 (ppb h)")) +
  stat_poly_eq(
    aes(label = paste(after_stat(rr.label),
      after_stat(n.label),
      sep = "*\", \"*"
    )),
    method = "lm", family = "Roboto Condensed", size = 5
  ) +
  # scale_x_continuous(limits = c(0, 23e3)) +
  # scale_y_continuous(limits = c(0, 23e3)) +
  geom_xsidehistogram(bins = 60) +
  geom_ysidehistogram(bins = 60) +
  ggside(x.pos = "top", y.pos = "right", collapse = "all") +
  scale_xsidey_continuous(labels = NULL, limits = c(0, 2e4)) +
  scale_ysidex_continuous(labels = NULL, limits = c(0, 2e4)) +
  theme_half_open(18, font_family = "Roboto Condensed") +
  background_grid() +
  theme(legend.position = c(0.6, 0.2))

patch <-
  # (p1 + p2) / (p3 + p4) +
  (p3 + p4) +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(size = 30, family = "Roboto Condensed"))

ggsave("figures_RB/daily_gam_performance.pdf", patch, width = 2, height = 1, scale = 8)
