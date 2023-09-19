source("script/loadPackages.R")
source("script/loadFunctions.R")
source("script/loadFormulas.R")

data <- get_data() %>%
  mutate(yeari = fcase(
    year %in% 2005:2009, "old",
    year %in% 2010:2014, "mid",
    year %in% 2015:2019, "new"
  ))

data <- data %>%
  nest(fdata = -crop_parent)

data$res <- map(data$fdata, function(adata) {
  adata <- adata %>%
    mutate(AOD2 = AOD^2)

  adata <- map(c(
    "AOD", "AOD2", "AOT40"
  ), function(avar) {
    i(adata$yeari, adata %>% pull(!!avar)) %>%
      as_tibble() %>%
      `names<-`(str_c(names(.), "_", avar))
  }) %>%
    bind_cols() %>%
    bind_cols(adata)

  model <- feols(fml_inter_year, adata,
    cluster = ~county,
    weights = ~fraction, nthreads = 0, mem.clean = T
  )

  AOD_range <- range(adata$AOD)
  AOD_hist <- seq(AOD_range[1], AOD_range[2], 0.01)
  AOD_mid <- round(quantile(adata$AOD, 0.5), 2)

  AOD_res <- map_dfr(c("old", "mid", "new"), function(atime) {
    AOD_tbl <- tibble(AOD_hist, AOD_hist^2) %>% `names<-`(str_c(atime, c("_AOD", "_AOD2")))
    AOD_tbl_mid <- tibble(rep(AOD_mid, length(AOD_hist)), rep(AOD_mid^2, length(AOD_hist))) %>% `names<-`(str_c(atime, c("_AOD", "_AOD2")))


    relpred(model, AOD_tbl, AOD_tbl_mid) %>%
      mutate(x = AOD_hist, yeari = atime)
  })

  ozone_range <- range(adata$AOT40)
  ozone_hist <- seq(ozone_range[1], ozone_range[2], 100)
  ozone_mid <- round(quantile(adata$AOT40, 0.5), 0)

  ozone_res <- map_dfr(c("old", "mid", "new"), function(atime) {
    ozone_tbl <- tibble(ozone_hist) %>% `names<-`(str_c(atime, c("_AOT40")))
    ozone_tbl_mid <- tibble(rep(ozone_mid, length(ozone_hist))) %>% `names<-`(str_c(atime, c("_AOT40")))


    relpred(model, ozone_tbl, ozone_tbl_mid) %>%
      mutate(x = ozone_hist, yeari = atime)
  })

  lst(AOD_res, ozone_res)
})

p1 <- data %>%
  select(-fdata) %>%
  unnest_wider(res) %>%
  select(crop_parent, ozone_res) %>%
  unnest() %>%
  mutate(
    yeari = factor(yeari, c("old", "mid", "new"), c("2005-2009", "2010-2014", "2015-2019")),
    across(c(fit, lwr, upr), ~ expm1(.x) * 1e2)
  ) %>%
  ggplot(aes(x = x, y = fit, color = yeari)) +
  facet_wrap(vars(crop_parent), scales = "free", ncol = 1) +
  geom_ribbon(aes(x = x, ymin = lwr, ymax = upr, fill = yeari),
    show.legend = F, alpha = 0.3, color = NA
  ) +
  geom_line(size = 2, show.legend = F) +
  geom_hline(yintercept = 0, size = 1, linetype = "dashed") +
  geom_xsidehistogram(
    aes(x = AOT40, y = ..ncount.., fill = yeari),
    data = data %>% select(-res) %>% unnest() %>% mutate(yeari = factor(yeari, c("old", "mid", "new"), c("2005-2009", "2010-2014", "2015-2019"))),
    show.legend = F, alpha = 0.6, color = "black",
    size = 0.3, bins = 60, position = "identity", inherit.aes = F
  ) +
  xlab("AOT40 (ppb h)") +
  theme_half_open(18, font_family = "Roboto Condensed") +
  background_grid() +
  theme(
    ggside.panel.scale = 0.2,
    legend.title = element_blank(),
    legend.text = element_text(size = 15)
  ) +
  ggside(x.pos = "bottom") +
  scale_x_continuous(expand = expansion()) +
  scale_xsidey_continuous(labels = NULL, breaks = NULL)

p2 <- data %>%
  select(-fdata) %>%
  unnest_wider(res) %>%
  select(crop_parent, AOD_res) %>%
  unnest() %>%
  mutate(
    yeari = factor(yeari, c("old", "mid", "new"), c("2005-2009", "2010-2014", "2015-2019")),
    across(c(fit, lwr, upr), ~ expm1(.x) * 1e2)
  ) %>%
  ggplot(aes(x = x, y = fit, color = yeari)) +
  facet_wrap(vars(crop_parent), scales = "free", ncol = 1) +
  geom_ribbon(aes(x = x, ymin = lwr, ymax = upr, fill = yeari),
    show.legend = F, alpha = 0.3, color = NA
  ) +
  geom_line(size = 2, show.legend = T) +
  geom_hline(yintercept = 0, size = 1, linetype = "dashed") +
  geom_xsidehistogram(
    aes(x = AOD, y = ..ncount.., fill = yeari),
    data = data %>% select(-res) %>% unnest() %>% mutate(yeari = factor(yeari, c("old", "mid", "new"), c("2005-2009", "2010-2014", "2015-2019"))),
    show.legend = F, alpha = 0.6, color = "black",
    size = 0.3, bins = 60, position = "identity", inherit.aes = F
  ) +
  xlab("AOD") +
  theme_half_open(18, font_family = "Roboto Condensed") +
  background_grid() +
  theme(
    ggside.panel.scale = 0.2,
    legend.title = element_blank(),
    legend.text = element_text(size = 15),
    legend.position = c(1.1, 0.5)
  ) +
  ggside(x.pos = "bottom") +
  scale_x_continuous(expand = expansion()) +
  scale_xsidey_continuous(labels = NULL, breaks = NULL)

patch <- p1 + p2 + plot_spacer() +
  plot_layout(nrow = 1, widths = c(1, 1, 0.5)) &
  ylab("Percentage Change in SIF") &
  scale_color_manual(values = c("#008500", "#D27685", "#FF8500")) &
  scale_fill_manual(values = c("#008500", "#D27685", "#FF8500"))

ggsave("figures/response_hete_temporal.pdf", patch, width = 2, height = 2.25, scale = 6)
