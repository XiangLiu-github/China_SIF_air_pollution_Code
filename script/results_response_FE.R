source("script/loadPackages.R")
source("script/loadFunctions.R")
source("script/loadFormulas.R")

data <- get_data()

data <- data %>%
  nest(fdata = -crop_parent)

model_names <- c(
  "preferred", "grid quadratic trend", "county year trend"
) %>%
  str_c("(", 1:length(.), ") ", .)

get_plot <- function(adata, acrop) {
  models <- lst(
    fml_base, fml_grid_year_trend2, fml_county_year_trend
  ) %>%
    map(~ feols(.x, adata,
      cluster = ~county,
      weights = ~fraction, nthreads = 0, lean = T
    ), .progress = T) %>%
    set_names(model_names)

  colors <- c("#FEE3EC", "#FFF9B6", "#EAF5FF")

  ozone_range <- range(adata$AOT40)
  ozone_hist <- seq(ozone_range[1], ozone_range[2], 100)
  ozone_mid <- round(quantile(adata$AOT40, 0.5), 0)
  ozone <- map_dfr(models, function(amodel) {
    relpred(
      amodel,
      tibble(AOT40 = ozone_hist),
      tibble(AOT40 = rep(ozone_mid, length(ozone_hist)))
    ) %>%
      mutate(ozone_x = ozone_hist)
  }, .id = "model") %>%
    mutate(across(c(fit, lwr, upr), ~ expm1(.x) * 1e2))

  ozone_offset <- (ozone_range[2] - ozone_range[1]) * 0.02

  ozone_lab <- map_dfr(models, function(amodel) {
    relpred(
      amodel,
      tibble(AOT40 = rep(ozone_range[2], 10)),
      tibble(AOT40 = rep(ozone_mid, length(rep(ozone_range[2], 10))))
    ) %>%
      mutate(ozone_x = rep(ozone_range[2], 10))
  }, .id = "model") %>%
    group_by(model) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    separate(model, c("label", NA), "\\s", remove = F) %>%
    mutate(across(c(fit, lwr, upr), ~ expm1(.x) * 1e2)) %>%
    mutate(ozone_x = ozone_x + ozone_offset)

  p1 <- ozone %>%
    mutate(model = fct_relevel(model, rev(model_names))) %>%
    ggplot(aes(x = ozone_x, y = fit, color = model)) +
    geom_hline(yintercept = 0, size = 1, linetype = "dashed") +
    geom_ribbon(aes(x = ozone_x, ymin = lwr, ymax = upr),
      inherit.aes = F,
      data = . %>% filter(str_detect(model, "preferred")),
      show.legend = F, fill = colors[1]
    ) +
    geom_line(size = 1, show.legend = F) +
    geom_text_repel(aes(label = label),
      data = ozone_lab,
      show.legend = F, min.segment.length = 0, box.padding = 1,
      direction = "y", xlim = c(ozone_range[2], NA), max.overlaps = Inf,
      arrow = arrow,
      family = "Roboto Condensed", size = 5, seed = 2022
    ) +
    # ylab(TeX("Change in SIF (mW $m^{-2}$ $nm^{-1}$ $sr^{-1}$)")) +
    ylab("Percentage Change in SIF") +
    xlab("AOT40 (ppb h)") +
    theme_half_open(18, font_family = "Roboto Condensed") +
    background_grid() +
    theme(
      ggside.panel.scale = 0.2,
      legend.title = element_blank(),
      legend.position = c(0.6, 0.8)
    ) +
    geom_xsidehistogram(aes(x = AOT40),
      data = adata,
      bins = 60, fill = colors[1], color = "black", size = 0.3,
      inherit.aes = F
    ) +
    ggside(x.pos = "bottom", collapse = "x") +
    scale_x_continuous(limits = c(ozone_range[1], ozone_range[2] + ozone_offset), expand = expansion(mult = c(0, 0.15))) +
    scale_xsidey_continuous(labels = NULL, breaks = NULL) +
    labs(subtitle = acrop)

  AOD_range <- range(adata$AOD)
  AOD_hist <- seq(AOD_range[1], AOD_range[2], 0.01)
  AOD_mid <- round(quantile(adata$AOD, 0.5), 2)
  AOD <- map_dfr(models, function(amodel) {
    relpred(
      amodel,
      tibble(AOD = AOD_hist),
      tibble(AOD = rep(AOD_mid, length(AOD_hist)))
    ) %>%
      mutate(AOD_x = AOD_hist)
  }, .id = "model") %>%
    mutate(across(c(fit, lwr, upr), ~ expm1(.x) * 1e2))

  AOD_offset <- (AOD_range[2] - AOD_range[1]) * 0.02

  AOD_lab <- map_dfr(models, function(amodel) {
    relpred(
      amodel,
      tibble(AOD = rep(AOD_range[2], 10)),
      tibble(AOD = rep(AOD_mid, length(rep(AOD_range[2], 10))))
    ) %>%
      mutate(AOD_x = rep(AOD_range[2], 10))
  }, .id = "model") %>%
    group_by(model) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    separate(model, c("label", NA), "\\s", remove = F) %>%
    mutate(across(c(fit, lwr, upr), ~ expm1(.x) * 1e2)) %>%
    mutate(AOD_x = AOD_x + AOD_offset)

  p2 <- AOD %>%
    mutate(model = fct_relevel(model, rev(model_names))) %>%
    ggplot(aes(x = AOD_x, y = fit, color = model)) +
    geom_hline(yintercept = 0, size = 1, linetype = "dashed") +
    # geom_vline(aes(xintercept = AOD_curving_point, ),
    #            data = NULL, size = 1, linetype = "dashed",
    #            show.legend = F
    # ) +
    geom_ribbon(aes(x = AOD_x, ymin = lwr, ymax = upr),
      inherit.aes = F,
      data = . %>% filter(str_detect(model, "preferred")),
      show.legend = F, fill = colors[2]
    ) +
    geom_line(size = 1, show.legend = T) +
    # geom_text(aes(x = x, y = y, label = label),
    #           data = tibble(
    #             x = AOD_curving_point - 0.04, y = min(AOD$fit),
    #             label = sprintf("%0.2f", round(AOD_curving_point, 2))
    #           ),
    #           size = 6,
    #           inherit.aes = F, show.legend = F, family = "Roboto Condensed"
    # ) +
    geom_text_repel(aes(label = label),
      data = AOD_lab,
      show.legend = F, min.segment.length = 0, box.padding = 1,
      direction = "y", xlim = c(AOD_range[2], NA), max.overlaps = Inf,
      arrow = arrow,
      family = "Roboto Condensed", size = 5, seed = 2022
    ) +
    # ylab(TeX("Change in SIF (mW $m^{-2}$ $nm^{-1}$ $sr^{-1}$)")) +
    ylab("Percentage Change in SIF") +
    xlab("AOD") +
    theme_half_open(18, font_family = "Roboto Condensed") +
    background_grid() +
    theme(
      ggside.panel.scale = 0.2,
      legend.title = element_blank(),
      legend.position = c(1.05, 0.5)
    ) +
    geom_xsidehistogram(aes(x = AOD),
      data = adata,
      bins = 60, fill = colors[2], color = "black", size = 0.3,
      inherit.aes = F
    ) +
    ggside(x.pos = "bottom", collapse = "x") +
    scale_x_continuous(limits = c(AOD_range[1], AOD_range[2] + AOD_offset), expand = expansion(mult = c(0, 0.15))) +
    scale_xsidey_continuous(labels = NULL, breaks = NULL) +
    labs(subtitle = acrop)

  return(lst(p1, p2, plot_spacer()))
}

data <- data %>%
  mutate(plots = map2(fdata, crop_parent, get_plot))

patch <- reduce(flatten(data$plots), `+`) +
  plot_layout(ncol = 3, widths = c(1, 1, 0.5)) &
  scale_color_manual(
    values = wesanderson::wes_palette("Rushmore1",
      n = length(model_names),
      type = c("continuous")
    ),
    guide = guide_legend(reverse = T)
  ) &
  theme(plot.subtitle = element_text(size = 20))

ggsave("figures/response_FE.pdf", patch, width = 3, height = 3, scale = 5)
