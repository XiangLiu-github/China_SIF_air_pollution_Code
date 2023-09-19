source("script/loadPackages.R")
source("script/loadFunctions.R")

data <- get_data(varss = c("GOSIF_sum", "down", "diff", "dir", "cloud", "AOD"))

models_pooled <- feols(
  c(down, diff, dir, dif) ~
    cloud + cloud^2 + cloud^3 +
      AOD + AOD^2 + AOD^3 |
      x_y^crop_parent + year^crop_parent,
  data = data, weights = ~fraction,
  cluster = ~county, nthreads = 0
) %>%
  as.list()

models_indi <- feols(
  c(down, diff, dir, dif) ~
    cloud + cloud^2 + cloud^3 +
      AOD + AOD^2 + AOD^3 |
      x_y + year,
  weights = ~fraction,
  data = data, split = ~crop_parent,
  cluster = ~county, nthreads = 0
) %>%
  as.list()

models <- append(models_pooled, models_indi)

crops <- c("Pooled", "Maize", "Rice", "Wheat")
rads <- c("Global radiation", "Diffuse radiation", "Direct radiation", "Diffuse fraction")

names(models) <-
  str_c(rep(crops, each = 4), "_", rep(rads, times = 4))

cloud_range <- range(data$cloud)
cloud_hist <- seq(cloud_range[1], cloud_range[2], 0.1)
cloud_mid <- round(quantile(data$cloud, 0.5), 0)

cloud_tbl <- map_dfr(models, function(amodel) {
  relpred(
    amodel,
    tibble(cloud = cloud_hist),
    tibble(cloud = rep(cloud_mid, length(cloud_hist)))
  ) %>%
    mutate(cloud_x = cloud_hist)
}, .id = "model") %>%
  separate(model, c("crop_parent", "model"), "_") %>%
  mutate(
    crop_parent = factor(crop_parent, crops),
    model = factor(model, c("Global radiation", "Direct radiation", "Diffuse radiation", "Diffuse fraction"))
  )

p1 <- cloud_tbl %>%
  ggplot(aes(x = cloud_x, y = fit, color = crop_parent)) +
  facet_wrap(~model, scales = "free_y", ncol = 1) +
  geom_ribbon(aes(x = cloud_x, ymin = lwr, ymax = upr, fill = crop_parent),
    show.legend = F, alpha = 0.3, color = NA
  ) +
  geom_line(size = 2, show.legend = F) +
  gghighlight(crop_parent == "Pooled",
    n = 1, use_direct_label = FALSE,
    keep_scales = TRUE, calculate_per_facet = TRUE,
    unhighlighted_params = list(colour = NULL, alpha = 0.3)
  ) +
  geom_hline(yintercept = 0, size = 1, linetype = "dashed") +
  xlab(latex2exp::TeX("COD")) +
  scale_x_continuous(limits = cloud_range, expand = expansion()) +
  theme_half_open(18, font_family = "Roboto Condensed") +
  background_grid()

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
  mutate(
    crop_parent = factor(crop_parent, crops),
    model = factor(model, c("Global radiation", "Direct radiation", "Diffuse radiation", "Diffuse fraction"))
  )

p2 <- AOD_tbl %>%
  ggplot(aes(x = AOD_x, y = fit, color = crop_parent)) +
  facet_wrap(~model, scales = "free_y", ncol = 1) +
  geom_ribbon(aes(x = AOD_x, ymin = lwr, ymax = upr, fill = crop_parent),
    show.legend = F, alpha = 0.3, color = NA
  ) +
  geom_line(size = 2, show.legend = T) +
  gghighlight(crop_parent == "Pooled",
    n = 1, use_direct_label = FALSE,
    keep_scales = TRUE, calculate_per_facet = TRUE,
    unhighlighted_params = list(colour = NULL, alpha = 0.3)
  ) +
  geom_hline(yintercept = 0, size = 1, linetype = "dashed") +
  xlab(latex2exp::TeX("AOD")) +
  scale_x_continuous(limits = AOD_range, expand = expansion()) +
  theme_half_open(18, font_family = "Roboto Condensed") +
  background_grid() +
  theme(
    legend.title = element_blank(),
    legend.position = c(1.1, 0.5)
  )

patch <- p1 + p2 + plot_spacer() +
  plot_layout(nrow = 1, widths = c(1, 1, 0.3)) &
  ylab(TeX("Change in radition components")) &
  scale_color_manual(values = NatParksPalettes::natparks.pals("Yellowstone", n = length(crops))) &
  scale_fill_manual(values = NatParksPalettes::natparks.pals("Yellowstone", n = length(crops)))

ggsave("figures/response_radiation_Chakraborty.pdf", patch, width = 2.3, height = 2, scale = 5)
