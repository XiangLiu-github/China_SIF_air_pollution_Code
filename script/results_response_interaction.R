source("script/loadPackages.R")
source("script/loadFunctions.R")
source("script/loadFormulas.R")

data <- get_data()

data <- data %>%
  group_by(crop, x, y) %>%
  mutate(vpd_mean = fmean(vpd_mean)) %>%
  group_by(crop_parent) %>%
  mutate(vpd_mean = DescTools::Winsorize(vpd_mean, minval = 3, maxval = 15)) %>%
  ungroup()

data <- data %>%
  nest(.by = crop_parent, .key = "fdata")

data <- data %>%
  mutate(
    model_v = map2(crop_parent, fdata, function(acrop, adata) {
      feols(fml_inter_v, adata,
        cluster = ~county,
        weights = ~fraction, nthreads = qn
      )
    }),
    model_i = map2(crop_parent, fdata, function(acrop, adata) {
      feols(fml_inter_i, adata,
        cluster = ~county,
        weights = ~fraction, nthreads = qn
      )
    })
  )

data_raw <- data %>%
  select(-starts_with("model")) %>%
  unnest()

ozone <- data %>%
  mutate(datas = pmap(list(fdata, model_v, model_i), function(afdata, amodel_v, amodel_i) {
    ozone_range <- round(range(afdata$AOT40))
    ozone_hist <- seq(ozone_range[1], ozone_range[2], length.out = 20)
    ozone_mid <- 1e4

    v_range <- round(range(afdata$vpd_mean), 1)
    v_hist <- seq(v_range[1], v_range[2], length.out = 20)

    v_res <- map_dfr(v_hist, function(av_hist) {
      v_crossing <- crossing(vpd_mean = av_hist, AOT40 = ozone_hist)
      v_crossing_mid <- v_crossing %>% mutate(AOT40 = ozone_mid)

      relpred(
        amodel_v,
        v_crossing,
        v_crossing_mid
      ) %>%
        bind_cols(v_crossing)
    })

    i_range <- c(0, 1)
    i_hist <- seq(i_range[1], i_range[2], length.out = 20)

    i_res <- map_dfr(i_hist, function(ai_hist) {
      i_crossing <- crossing(irg_fraction = ai_hist, AOT40 = ozone_hist)
      i_crossing_mid <- i_crossing %>% mutate(AOT40 = ozone_mid)

      relpred(
        amodel_i,
        i_crossing,
        i_crossing_mid
      ) %>%
        bind_cols(i_crossing)
    })

    lst(v_res, i_res)
  }), .keep = "unused")

p1 <- ozone %>%
  unnest_wider(datas) %>%
  select(-i_res) %>%
  unnest() %>%
  mutate(across(c(fit, lwr, upr), ~ expm1(.x) * 1e2)) %>%
  ggplot(aes(x = AOT40, y = vpd_mean)) +
  facet_wrap(~crop_parent, scales = "free") +
  geom_raster(aes(fill = fit),
    hjust = 0,
    vjust = 0,
    interpolate = T
  ) +
  geom_point(data = . %>% filter(upr <= 0 | lwr >= 0), color = "grey", size = 1, shape = 15) +
  geom_vline(xintercept = 1e4, linetype = "longdash") +
  geom_xsidehistogram(aes(x = AOT40, y = ..ncount..),
    data = data_raw,
    bins = 60, color = "black", size = 0.3,
    inherit.aes = F
  ) +
  geom_ysidehistogram(aes(y = vpd_mean, x = ..ncount..),
    data = data_raw,
    bins = 60, color = "black", size = 0.3,
    inherit.aes = F
  ) +
  scale_y_continuous(name = "VPD (hPa)", expand = expansion()) +
  scale_x_continuous(name = "AOT40 (ppb h)", expand = expansion())

p2 <- ozone %>%
  filter(crop_parent != "Rice") %>%
  unnest_wider(datas) %>%
  select(-v_res) %>%
  unnest() %>%
  mutate(across(c(fit, lwr, upr), ~ expm1(.x) * 1e2)) %>%
  ggplot(aes(x = AOT40, y = irg_fraction)) +
  facet_wrap(~crop_parent, scales = "free") +
  geom_raster(aes(fill = fit),
    hjust = 0,
    vjust = 0,
    interpolate = T
  ) +
  geom_point(data = . %>% filter(upr <= 0 | lwr >= 0), color = "grey", size = 1, shape = 15) +
  geom_vline(xintercept = 1e4, linetype = "longdash") +
  geom_xsidehistogram(aes(x = AOT40, y = ..ncount..),
    data = data_raw %>% filter(crop_parent != "Rice"),
    bins = 60, color = "black", size = 0.3,
    inherit.aes = F
  ) +
  geom_ysidehistogram(aes(y = irg_fraction, x = ..ncount..),
    data = data_raw %>% filter(crop_parent != "Rice"),
    bins = 60, color = "black", size = 0.3,
    inherit.aes = F
  ) +
  scale_y_continuous(name = "Irrigation fraction", expand = expansion()) +
  scale_x_continuous(name = "AOT40 (ppb h)", expand = expansion())


AOD <- data %>%
  mutate(datas = pmap(list(fdata, model_v, model_i), function(afdata, amodel_v, amodel_i) {
    AOD_range <- round(range(afdata$AOD), 2)
    AOD_hist <- seq(AOD_range[1], AOD_range[2], length.out = 20)
    AOD_mid <- 0.5

    v_range <- round(range(afdata$vpd_mean), 1)
    v_hist <- seq(v_range[1], v_range[2], length.out = 20)

    v_res <- map_dfr(v_hist, function(av_hist) {
      v_crossing <- crossing(vpd_mean = av_hist, AOD = AOD_hist)
      v_crossing_mid <- v_crossing %>% mutate(AOD = AOD_mid)

      relpred(
        amodel_v,
        v_crossing,
        v_crossing_mid
      ) %>%
        bind_cols(v_crossing)
    })

    i_range <- c(0, 1)
    i_hist <- seq(i_range[1], i_range[2], length.out = 20)

    i_res <- map_dfr(i_hist, function(ai_hist) {
      i_crossing <- crossing(irg_fraction = ai_hist, AOD = AOD_hist)
      i_crossing_mid <- i_crossing %>% mutate(AOD = AOD_mid)

      relpred(
        amodel_i,
        i_crossing,
        i_crossing_mid
      ) %>%
        bind_cols(i_crossing)
    })

    lst(v_res, i_res)
  }), .keep = "unused")

p3 <- AOD %>%
  unnest_wider(datas) %>%
  select(-i_res) %>%
  unnest() %>%
  mutate(across(c(fit, lwr, upr), ~ expm1(.x) * 1e2)) %>%
  ggplot(aes(x = AOD, y = vpd_mean)) +
  facet_wrap(~crop_parent, scales = "free") +
  geom_raster(aes(fill = fit),
    hjust = 0,
    vjust = 0,
    interpolate = T
  ) +
  geom_point(data = . %>% filter(upr <= 0 | lwr >= 0), color = "grey", size = 1, shape = 15) +
  geom_vline(xintercept = 0.5, linetype = "longdash") +
  geom_xsidehistogram(aes(x = AOD, y = ..ncount..),
    data = data_raw,
    bins = 60, color = "black", size = 0.3,
    inherit.aes = F
  ) +
  geom_ysidehistogram(aes(y = vpd_mean, x = ..ncount..),
    data = data_raw,
    bins = 60, color = "black", size = 0.3,
    inherit.aes = F
  ) +
  scale_y_continuous(name = "VPD (hPa)", expand = expansion()) +
  scale_x_continuous(name = "AOD", expand = expansion())


p4 <- AOD %>%
  filter(crop_parent != "Rice") %>%
  unnest_wider(datas) %>%
  select(-v_res) %>%
  unnest() %>%
  mutate(across(c(fit, lwr, upr), ~ expm1(.x) * 1e2)) %>%
  ggplot(aes(x = AOD, y = irg_fraction)) +
  facet_wrap(~crop_parent, scales = "free") +
  geom_raster(aes(fill = fit),
    hjust = 0,
    vjust = 0,
    interpolate = T
  ) +
  geom_point(data = . %>% filter(upr <= 0 | lwr >= 0), color = "grey", size = 1, shape = 15) +
  geom_vline(xintercept = 0.5, linetype = "longdash") +
  geom_xsidehistogram(aes(x = AOD, y = ..ncount..),
    data = data_raw %>% filter(crop_parent != "Rice"),
    bins = 60, color = "black", size = 0.3,
    inherit.aes = F
  ) +
  geom_ysidehistogram(aes(y = irg_fraction, x = ..ncount..),
    data = data_raw %>% filter(crop_parent != "Rice"),
    bins = 60, color = "black", size = 0.3,
    inherit.aes = F
  ) +
  scale_y_continuous(name = "Irrigation fraction", expand = expansion()) +
  scale_x_continuous(name = "AOD", expand = expansion())


make_plot <- function(ap) {
  ap +
    scale_xsidey_continuous(labels = NULL, breaks = NULL) +
    scale_ysidex_continuous(labels = NULL, breaks = NULL) +
    theme_half_open(18, font_family = "Roboto Condensed") +
    background_grid("none") +
    theme(
      ggside.panel.scale = 0.2,
      legend.text = element_text(size = 15),
      legend.justification = "center",
      legend.key.height = unit(3, "lines"),
      legend.title = element_text(angle = -90, family = "Roboto Condensed")
    ) +
    ggside(x.pos = "bottom", y.pos = "left") +
    scale_fill_gradient2(
      name = "Percentage Change in SIF",
      limits = c(-10, 10),
      oob = squish
    ) +
    guides(fill = guide_colorbar(title.position = "right"))
}

patch1 <- list(p1, p3) %>%
  map(make_plot) %>%
  wrap_plots(ncol = 1) +
  plot_layout(guides = "collect")

ggsave("figures/response_hete_vpd.pdf", patch1, width = 3.5, height = 2.5, scale = 4)

patch2 <- list(p2, p4) %>%
  map(make_plot) %>%
  wrap_plots(ncol = 1) +
  plot_layout(guides = "collect")

ggsave("figures/response_hete_irrigation.pdf", patch2, width = 2.5, height = 2.5, scale = 4)
