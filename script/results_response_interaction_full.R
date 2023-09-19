source("script/loadPackages.R")
source("script/loadFunctions.R")
source("script/loadFormulas.R")

data <- get_data()

data <- data %>%
  group_by(crop, x, y) %>%
  mutate(vpd_mean = fmean(vpd_mean)) %>%
  group_by(crop_parent) %>%
  mutate(vpd_mean = DescTools::Winsorize(vpd_mean, probs = c(0.01, 0.99))) %>%
  ungroup()

data <- data %>%
  nest(fdata = -crop_parent)

data$res <- map2(data$crop_parent, data$fdata, function(acrop, adata) {
  adata <- adata %>%
    mutate(AOD2 = AOD^2)

  adata <- map(c(
    "AOD", "AOD2"
  ), function(avar) {
    i(adata$region, adata %>% pull(!!avar)) %>%
      as_tibble() %>%
      `names<-`(str_c(names(.), "_", avar))
  }) %>%
    bind_cols() %>%
    bind_cols(adata)

  if (acrop == "Rice") {
    model <- feols(fml_inter_full_rice, adata,
      cluster = ~county,
      weights = ~fraction, nthreads = 0, mem.clean = T
    )
  } else {
    model <- feols(fml_inter_full, adata,
      cluster = ~county,
      weights = ~fraction, nthreads = 0, mem.clean = T
    )
  }

  AOD_range <- range(adata$AOD)
  AOD_hist <- seq(AOD_range[1], AOD_range[2], 0.01)
  AOD_mid <- round(quantile(adata$AOD, 0.5), 2)

  AOD_res <- map_dfr(unique(adata$region), function(atime) {
    AOD_tbl <- tibble(AOD_hist, AOD_hist^2) %>% `names<-`(str_c(atime, c("_AOD", "_AOD2")))
    AOD_tbl_mid <- tibble(rep(AOD_mid, length(AOD_hist)), rep(AOD_mid^2, length(AOD_hist))) %>% `names<-`(str_c(atime, c("_AOD", "_AOD2")))


    relpred(model, AOD_tbl, AOD_tbl_mid) %>%
      mutate(x = AOD_hist, region = atime)
  })

  ozone_range <- round(range(adata$AOT40))
  ozone_hist <- seq(ozone_range[1], ozone_range[2], length.out = 20)
  ozone_mid <- 1e4

  v_range <- round(range(adata$vpd_mean), 1)
  v_hist <- seq(v_range[1], v_range[2], length.out = 20)

  v_res <- map_dfr(v_hist, function(av_hist) {
    v_crossing <- crossing(vpd_mean = av_hist, AOT40 = ozone_hist)
    v_crossing_mid <- v_crossing %>% mutate(AOT40 = ozone_mid)

    relpred(
      model,
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
      model,
      i_crossing,
      i_crossing_mid
    ) %>%
      bind_cols(i_crossing)
  })

  ozone_res <- lst(v_res, i_res)

  lst(AOD_res, ozone_res)
})

data_raw <- data %>%
  select(-starts_with("res")) %>%
  unnest()

p1 <- data %>%
  unnest_wider(res) %>%
  select(crop_parent, ozone_res) %>%
  unnest_wider(ozone_res) %>%
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

p2 <- data %>%
  unnest_wider(res) %>%
  select(crop_parent, ozone_res) %>%
  unnest_wider(ozone_res) %>%
  select(-v_res) %>%
  filter(crop_parent != "Rice") %>%
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
      limits = c(-20, 20),
      oob = squish
    ) +
    guides(fill = guide_colorbar(title.position = "right"))
}

patch1 <- list(p1, p2) %>%
  map(make_plot) %>%
  wrap_plots(ncol = 1) +
  plot_layout(guides = "collect")
