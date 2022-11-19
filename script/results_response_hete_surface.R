source("script/loadPackages.R")
source("script/loadFunctions.R")
source("script/loadFormulas.R")

data <- get_data()

data <- data %>%
  group_by(x_y) %>%
  mutate(surface_mean = mean(surface)) %>%
  group_by(crop_parent) %>%
  mutate(surface_mean = fcase(
    surface_mean >= quantile(surface, 0.75), "wet",
    surface_mean <= quantile(surface, 0.25), "dry"
  )) %>%
  ungroup() %>%
  drop_na(surface_mean) %>%
  mutate(class = str_c(crop_parent, "_", surface_mean))

models <-
  feols(fml_base,
    data = data, split = ~class, weights = ~fraction,
    cluster = ~county, nthreads = 0, lean = T
  )

ozone_range <- range(data$AOT40)
ozone_hist <- seq(ozone_range[1], ozone_range[2], 100)
ozone_mid <- round(quantile(data$AOT40, 0.5), 0)
ozone <- map_dfr(models, function(amodel) {
  relpred(
    amodel,
    tibble(AOT40 = ozone_hist),
    tibble(AOT40 = rep(ozone_mid, length(ozone_hist)))
  ) %>%
    mutate(ozone_x = ozone_hist)
}, .id = "model") %>%
  separate(model, c("crop_parent", "model"), "_") %>% 
  mutate(across(c(fit, lwr, upr), ~ expm1(.x) * 1e2))

p1 <- ozone %>%
  ggplot(aes(x = ozone_x, y = fit, color = model)) +
  facet_wrap(vars(crop_parent), scales = "free", ncol = 1) +
  geom_ribbon(aes(x = ozone_x, ymin = lwr, ymax = upr, fill = model),
    show.legend = F, alpha = 0.3, color = NA
  ) +
  geom_line(size = 2, show.legend = F) +
  geom_hline(yintercept = 0, size = 1, linetype = "dashed") +
  geom_xsidehistogram(
    aes(x = AOT40, y = ..ncount.., fill = surface_mean),
    data = data,
    show.legend = F, alpha = 0.6, color = "black",
    size = 0.3, bins = 60, position = "identity", inherit.aes = F
  ) +
  xlab("AOT40 (ppb h)") +
  scale_x_continuous(limits = ozone_range, expand = expansion()) +
  theme_half_open(18, font_family = "Roboto Condensed") +
  background_grid() +
  theme(
    ggside.panel.scale = 0.2,
    legend.title = element_blank(),
    legend.text = element_text(size = 15)
  ) +
  ggside(x.pos = "bottom")

AOD_range <- range(data$AOD)
AOD_hist <- seq(AOD_range[1], AOD_range[2], 0.01)
AOD_mid <- round(quantile(data$AOD, 0.5), 2)
AOD_tbl <- map_dfr(models, function(amodel) {
  relpred(
    amodel,
    tibble(AOD = AOD_hist),
    tibble(AOD = rep(AOD_mid, length(AOD_hist)))
  ) %>%
    mutate(AOD_x = AOD_hist)
}, .id = "model") %>%
  separate(model, c("crop_parent", "model"), "_") %>% 
  mutate(across(c(fit, lwr, upr), ~ expm1(.x) * 1e2))

p2 <- AOD_tbl %>%
  ggplot(aes(x = AOD_x, y = fit, color = model)) +
  facet_wrap(vars(crop_parent), scales = "free", ncol = 1) +
  geom_ribbon(aes(x = AOD_x, ymin = lwr, ymax = upr, fill = model),
    show.legend = F, alpha = 0.3, color = NA
  ) +
  geom_line(size = 2, show.legend = T) +
  geom_hline(yintercept = 0, size = 1, linetype = "dashed") +
  geom_xsidehistogram(
    aes(x = AOD, y = ..ncount.., fill = surface_mean),
    data = data,
    show.legend = F, alpha = 0.6, color = "black",
    size = 0.3, bins = 60, position = "identity", inherit.aes = F
  ) +
  xlab(latex2exp::TeX("AOD")) +
  scale_x_continuous(limits = AOD_range, expand = expansion()) +
  theme_half_open(18, font_family = "Roboto Condensed") +
  background_grid() +
  theme(
    ggside.panel.scale = 0.2,
    legend.title = element_blank(),
    legend.text = element_text(size = 15),
    legend.position = c(1.1, 0.5)
  ) +
  ggside(x.pos = "bottom")

patch <- p1 + p2 + plot_spacer() +
  plot_layout(nrow = 1, widths = c(1, 1, 0.3)) &
  # ylab(TeX("Change in SIF (mW $m^{-2}$ $nm^{-1}$ $sr^{-1}$)")) &
  ylab("Percentage Change in SIF") &
  scale_color_manual(values = c("#126883", "#9FA3FE")) &
  scale_fill_manual(values = c("#126883", "#9FA3FE"))

ggsave("figures/response_hete_surface.pdf", patch, width = 2, height = 2.25, scale = 6)
