source("script/loadPackages.R")
source("script/loadFunctions.R")

require(parallel)
nc <- 9 ## cluster size, set for example portability
if (detectCores() > 1) { ## no point otherwise
  cl <- makeCluster(nc)
  ## could also use makeForkCluster, but read warnings first!
} else {
  cl <- NULL
}

data <-
  inner_join(
    read_rds("../Vegetation_ozone/data/inputs/station_raw/ozone/China.rds"),
    read_rds("../pre_NS2021/data/inputs/station_raw/ozone/China.rds")
  ) %>%
  mutate(O3 = rollingo3 * 1e3 * 2)

model_w126 <- bam(W126 ~ s(O3), data = data, family = tw(), chunk.size = 5000, cluster = cl)
model_aot40 <- bam(AOT40 ~ s(O3), data = data, family = tw(), chunk.size = 5000, cluster = cl)

folds <- vfold_cv(data)

cv_w126 <- map_dfr(folds$splits, function(asplit) {
  model <- bam(W126 ~ s(O3),
    data = asplit %>% analysis(),
    family = tw(), chunk.size = 5000, cluster = cl
  )

  tibble(
    predicted = exp(predict(model, asplit %>% assessment(), cluster = cl)),
    raw = asplit %>% assessment() %>% pull(W126)
  )
}, .progress = T)

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

p3 <- data %>%
  ggplot(aes(x = O3, y = W126)) +
  geom_pointdensity(adjust = 4, method = "default") +
  scale_color_viridis(
    trans = "log10",
    name = "n_neighbors",
    limits = c(1, 4e2), oob = scales::squish
  ) +
  xlab(TeX("MDA8 ($\\mu g $ $ m ^{-3}$)")) +
  ylab(TeX("W126 (ppb h)")) +
  scale_y_continuous(limits = c(0, 35e3)) +
  scale_x_continuous(limits = c(0, 230)) +
  geom_xsidehistogram(bins = 60) +
  geom_ysidehistogram(bins = 60) +
  ggside(x.pos = "top", y.pos = "right", collapse = "all") +
  scale_xsidey_continuous(labels = NULL) +
  scale_ysidex_continuous(labels = NULL) +
  geom_line(aes(x = O3, y = W126_pred),
    inherit.aes = F,
    data = tibble(O3 = seq(0, 230, by = 0.1)) %>%
      mutate(W126_pred = exp(predict(model_w126, .))),
    size = 1
  ) +
  theme_half_open(18, font_family = "Roboto Condensed") +
  background_grid() +
  theme(legend.position = c(0.2, 0.6))

p4 <- cv_w126 %>%
  ggplot(aes(x = raw, y = predicted)) +
  geom_pointdensity(adjust = 4, method = "default") +
  scale_color_viridis(
    trans = "log10",
    name = "n_neighbors",
    limits = c(1, 1e2), oob = scales::squish
  ) +
  geom_abline(slope = 1, intercept = 0, size = 1, color = "blue") +
  geom_smooth(method = "lm", se = F, color = "red", size = 1) +
  xlab(TeX("Observation W126 (ppb h)")) +
  ylab(TeX("Prediction W126 (ppb h)")) +
  stat_poly_eq(
    aes(label = paste(after_stat(rr.label),
      after_stat(n.label),
      sep = "*\", \"*"
    )),
    method = "lm", family = "Roboto Condensed", size = 5
  ) +
  scale_x_continuous(limits = c(0, 35e3)) +
  scale_y_continuous(limits = c(0, 35e3)) +
  geom_xsidehistogram(bins = 60) +
  geom_ysidehistogram(bins = 60) +
  ggside(x.pos = "top", y.pos = "right", collapse = "all") +
  scale_xsidey_continuous(labels = NULL, limits = c(0, 2e4)) +
  scale_ysidex_continuous(labels = NULL, limits = c(0, 2e4)) +
  theme_half_open(18, font_family = "Roboto Condensed") +
  background_grid() +
  theme(legend.position = c(0.6, 0.2))

p1 <- data %>%
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
  scale_y_continuous(limits = c(0, 23e3)) +
  scale_x_continuous(limits = c(0, 230)) +
  geom_xsidehistogram(bins = 60) +
  geom_ysidehistogram(bins = 60) +
  ggside(x.pos = "top", y.pos = "right", collapse = "all") +
  scale_xsidey_continuous(labels = NULL) +
  scale_ysidex_continuous(labels = NULL) +
  geom_line(aes(x = O3, y = AOT40_pred),
    inherit.aes = F,
    data = tibble(O3 = seq(0, 230, by = 0.1)) %>%
      mutate(AOT40_pred = exp(predict(model_aot40, .))),
    size = 1
  ) +
  theme_half_open(18, font_family = "Roboto Condensed") +
  background_grid() +
  theme(legend.position = c(0.2, 0.6))

p2 <- cv_aot40 %>%
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
  scale_x_continuous(limits = c(0, 23e3)) +
  scale_y_continuous(limits = c(0, 23e3)) +
  geom_xsidehistogram(bins = 60) +
  geom_ysidehistogram(bins = 60) +
  ggside(x.pos = "top", y.pos = "right", collapse = "all") +
  scale_xsidey_continuous(labels = NULL, limits = c(0, 2e4)) +
  scale_ysidex_continuous(labels = NULL, limits = c(0, 2e4)) +
  theme_half_open(18, font_family = "Roboto Condensed") +
  background_grid() +
  theme(legend.position = c(0.6, 0.2))

patch <-
  (p1 + p2) / (p3 + p4) +
    # (p3 + p4) +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(size = 30, family = "Roboto Condensed"))

ggsave("figures/gam_performance.pdf", patch, width = 2, height = 2, scale = 8)
