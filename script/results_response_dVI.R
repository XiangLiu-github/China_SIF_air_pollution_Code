source("script/loadPackages.R")
source("script/loadFunctions.R")

data <- get_data() %>%
  rename_with(~ str_c(.x, "_mean"), c(NIRv, NDVI, EVI))

fml <- str_c("c(NIRv_mean, NDVI_mean, EVI_mean,
               NIRv_peak, NDVI_peak, EVI_peak,
               NIRv_sum, NDVI_sum, EVI_sum) ~
            AOD + AOD ^ 2 +
            AOT40 +
            cloud + cloud ^ 2 +
            ", add_terms(c(
  str_c("bin", 0:39),
  str_c("surface_", 1:9)
)), "| x_y[year]") %>%
  as.formula()

models <- feols(fml,
  data = data, split = ~crop_parent, weights = ~fraction,
  cluster = ~county, nthreads = 0, lean = T
)


# hist(data %>% filter(crop_parent == 'Rice') %>% pull(AOD))

names(models) <- names(models) %>%
  str_remove("sample.var: crop_parent; sample: ") %>%
  str_remove(" lhs: ")

crops <- sort(unique(data$crop_parent))

ozone_range <- range(data$AOT40)
ozone_hist <- seq(ozone_range[1], ozone_range[2], 100)
ozone_mid <- round(quantile(data$AOT40, 0.5), 0)
ozone_tbl <- map_dfr(models, function(amodel) {
  relpred(
    amodel,
    tibble(AOT40 = ozone_hist),
    tibble(AOT40 = rep(ozone_mid, length(ozone_hist)))
  ) %>%
    mutate(ozone_x = ozone_hist)
}, .id = "model") %>%
  separate(model, c("crop_parent", "model"), ";") %>%
  separate(model, c("model", "arm"), "_") %>%
  mutate(
    crop_parent = factor(crop_parent, crops),
    model = factor(model, c("NIRv", "NDVI", "EVI")),
    arm = factor(arm, c("mean", "peak", "sum"))
  )

p1 <- ozone_tbl %>%
  filter(model != "Wenetal") %>%
  ggplot(aes(x = ozone_x, y = fit, color = model)) +
  facet_grid2(vars(arm), vars(crop_parent), scales = "free", independent = "y") +
  geom_ribbon(aes(x = ozone_x, ymin = lwr, ymax = upr, fill = model),
    show.legend = F, alpha = 0.3, color = NA
  ) +
  geom_line(size = 2, show.legend = F) +
  geom_hline(yintercept = 0, size = 1, linetype = "dashed") +
  xlab("AOT40 (ppb h)") +
  scale_x_continuous(limits = ozone_range, expand = expansion()) +
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
  separate(model, c("crop_parent", "model"), ";") %>%
  separate(model, c("model", "arm"), "_") %>%
  mutate(
    crop_parent = factor(crop_parent, crops),
    model = factor(model, c("NIRv", "NDVI", "EVI")),
    arm = factor(arm, c("mean", "peak", "sum"))
  )

p2 <- AOD_tbl %>%
  filter(model != "Wenetal") %>%
  ggplot(aes(x = AOD_x, y = fit, color = model)) +
  facet_grid2(vars(arm), vars(crop_parent), scales = "free", independent = "y") +
  geom_ribbon(aes(x = AOD_x, ymin = lwr, ymax = upr, fill = model),
    show.legend = F, alpha = 0.3, color = NA
  ) +
  geom_line(size = 2, show.legend = T) +
  geom_hline(yintercept = 0, size = 1, linetype = "dashed") +
  xlab(latex2exp::TeX("AOD")) +
  scale_x_continuous(limits = AOD_range, expand = expansion()) +
  theme_half_open(18, font_family = "Roboto Condensed") +
  background_grid() +
  theme(
    legend.title = element_blank(),
    legend.position = c(1.1, 1.3)
  )

patch <- p1 + p2 + plot_spacer() +
  plot_layout(
    design =
      "AAAAAC
                 BBBBBC"
  ) &
  ylab("Change in vegetation index") &
  scale_color_manual(values = c("#F17256", "#FCEFB2", "#3AC2DD")) &
  scale_fill_manual(values = c("#F17256", "#FCEFB2", "#3AC2DD"))

ggsave("figures/response_dVI.pdf", patch, width = 1, height = 1, scale = 13)
