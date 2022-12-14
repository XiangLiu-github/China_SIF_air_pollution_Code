source("script/loadPackages.R")
source("script/loadFunctions.R")
source("script/loadFormulas.R")

varss =
  c(
    "GOSIF_sum",
    "cloud",
    "AOD",
    "PM25",
    "W126_1",
    str_c("bin", 0:39),
    str_c("step", 1:10),
    str_c("root_", 1:9)
  )

data <- qread("data/tidied.qs", nthreads = qn) %>%
  lazy_dt() %>%
  # filter(GOSIF_sum > 0.05 & CSIF_sum > 0.05 & RTSIF_sum > 0.05) %>%
  mutate(
    crop_parent = fifelse(str_detect(crop, "Rice"), "Rice", crop),
    x_y = str_c(crop, x_y),
    across(
      c(
        starts_with("GOSIF"),
        starts_with("Wenetal"),
        starts_with("RTSIF"),
        starts_with("CSIF")
      ),
      log
    ),
    DF = DR / GR,
    DirR = GR - DR,
    AOT40 = AOT40_7,
    W126 = W126_7,
    AOT40 = fcase(
      crop_parent == "Maize",
      AOT40_7,
      crop_parent == "Rice",
      AOT40_7,
      crop_parent == "Wheat",
      AOT40_7 - AOT40_1
    ),
    W126 = fcase(
      crop_parent == "Maize",
      W126_7,
      crop_parent == "Rice",
      W126_7,
      crop_parent == "Wheat",
      W126_7 - W126_1
    ),
    O3 = fcase(
      crop_parent == "Maize",
      O3_7,
      crop_parent == "Rice",
      O3_7,
      crop_parent == "Wheat",
      (O3_7 * (MA - `GR&EM` + 1) - O3_1) / (MA - `GR&EM` + 0)
    )
    # O3 = O3_7
  ) %>%
  drop_na(all_of(varss)) %>%
  filter(province %in% c("河南省", "山东省", "山西省", "北京市", "天津市",
                         "安徽省", "河北省", "江苏省"),
         crop_parent != "Rice") %>% 
  group_by(crop, x, y) %>%
  filter(n() >= 10) %>%
  ungroup() %>%
  group_by(crop) %>%
  mutate(across(
    c(AOD, starts_with("W126"), starts_with("AOT40")),
    ~ DescTools::Winsorize(.x, probs = c(0, 0.995), na.rm = T)
  )) %>%
  ungroup() %>%
  arrange(crop) %>%
  as_tibble()

data <- data %>%
  nest(fdata = -crop_parent)

model_names <- c(
  "preferred", "linear", "cubic polynomial",
  "cubic spline (3 degree)", "cubic spline (5 degree)",
  "no aerosol", "no ozone"
) %>%
  str_c("(", 1:length(.), ") ", .)

get_plot <- function(adata, acrop) {
  PM25_ns3 <- ns(adata$PM25, df = 3)
  AOT40_ns3 <- ns(adata$AOT40, df = 3)
  AOD_ns3 <- ns(adata$AOD, df = 3)
  
  PM25_ns5 <- ns(adata$PM25, df = 5)
  AOT40_ns5 <- ns(adata$AOT40, df = 5)
  AOD_ns5 <- ns(adata$AOD, df = 5)
  
  data_modeling <- bind_cols(
    adata,
    as_tibble(PM25_ns3) %>%
      set_names(~ str_c("PM25_ns3_", .x)),
    as_tibble(AOT40_ns3) %>%
      set_names(~ str_c("AOT40_ns3_", .x)),
    as_tibble(AOD_ns3) %>%
      set_names(~ str_c("AOD_ns3_", .x)),
    as_tibble(PM25_ns5) %>%
      set_names(~ str_c("PM25_ns5_", .x)),
    as_tibble(AOT40_ns5) %>%
      set_names(~ str_c("AOT40_ns5_", .x)),
    as_tibble(AOD_ns5) %>%
      set_names(~ str_c("AOD_ns5_", .x))
  )
  
  models <- lst(
    fml_base, fml_linear, fml_poly3, fml_ns3, fml_ns5,
    fml_noaerosol, fml_noozone
  ) %>%
    pro_map(~ feols(.x, data_modeling,
                    cluster = ~county,
                    weights = ~fraction, nthreads = 0, mem.clean = T
    )) %>%
    set_names(model_names)
  
  colors <- c("#FEE3EC", "#FFF9B6", "#EAF5FF")
  
  ozone_range <- range(data_modeling$AOT40)
  ozone_hist <- seq(ozone_range[1], ozone_range[2], 100)
  ozone_mid <- round(quantile(data_modeling$AOT40, 0.5), 0)
  ozone <- map_dfr(models, function(amodel) {
    relpred(
      amodel,
      tibble(AOT40 = ozone_hist) %>%
        bind_cols(
          as_tibble(predict(AOT40_ns3, ozone_hist)) %>%
            set_names(~ str_c("AOT40_ns3_", .x)),
          as_tibble(predict(AOT40_ns5, ozone_hist)) %>%
            set_names(~ str_c("AOT40_ns5_", .x))
        ),
      tibble(AOT40 = rep(ozone_mid, length(ozone_hist))) %>%
        bind_cols(
          as_tibble(predict(
            AOT40_ns3,
            rep(ozone_mid, length(ozone_hist))
          )) %>%
            set_names(~ str_c("AOT40_ns3_", .x)),
          as_tibble(predict(
            AOT40_ns5,
            rep(ozone_mid, length(ozone_hist))
          )) %>%
            set_names(~ str_c("AOT40_ns5_", .x))
        )
    ) %>%
      mutate(ozone_x = ozone_hist)
  }, .id = "model") %>% 
    mutate(across(c(fit, lwr, upr), ~ expm1(.x) * 1e2))
  
  ozone_lab <- map_dfr(models, function(amodel) {
    relpred(
      amodel,
      tibble(AOT40 = rep(ozone_range[2], 10)) %>%
        bind_cols(
          as_tibble(predict(AOT40_ns3, rep(ozone_range[2], 10))) %>%
            set_names(~ str_c("AOT40_ns3_", .x)),
          as_tibble(predict(AOT40_ns5, rep(ozone_range[2], 10))) %>%
            set_names(~ str_c("AOT40_ns5_", .x))
        ),
      tibble(AOT40 = rep(ozone_mid, length(rep(ozone_range[2], 10)))) %>%
        bind_cols(
          as_tibble(predict(
            AOT40_ns3,
            rep(ozone_mid, length(rep(ozone_range[2], 10)))
          )) %>%
            set_names(~ str_c("AOT40_ns3_", .x)),
          as_tibble(predict(
            AOT40_ns5,
            rep(ozone_mid, length(rep(ozone_range[2], 10)))
          )) %>%
            set_names(~ str_c("AOT40_ns5_", .x))
        )
    ) %>%
      mutate(ozone_x = rep(ozone_range[2], 10))
  }, .id = "model") %>%
    group_by(model) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    separate(model, c("label", NA), "\\s", remove = F)  %>% 
    mutate(across(c(fit, lwr, upr), ~ expm1(.x) * 1e2))
  
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
                        data = data_modeling,
                        bins = 60, fill = colors[1], color = "black", size = 0.3,
                        inherit.aes = F
    ) +
    ggside(x.pos = "bottom", collapse = "x") +
    scale_x_continuous(limits = ozone_range, expand = expansion(mult = c(0, 0.15))) +
    scale_xsidey_continuous(labels = NULL) +
    labs(subtitle = acrop)
  
  # AOD_curving_point <- -coef(models$preferred)["AOD"] / (2 * coef(models$preferred)["I(AOD^2)"])
  
  AOD_range <- range(data_modeling$AOD)
  AOD_hist <- seq(AOD_range[1], AOD_range[2], 0.01)
  AOD_mid <- round(quantile(data_modeling$AOD, 0.5), 2)
  AOD <- map_dfr(models, function(amodel) {
    relpred(
      amodel,
      tibble(AOD = AOD_hist) %>%
        bind_cols(
          as_tibble(predict(AOD_ns3, AOD_hist)) %>%
            set_names(~ str_c("AOD_ns3_", .x)),
          as_tibble(predict(AOD_ns5, AOD_hist)) %>%
            set_names(~ str_c("AOD_ns5_", .x))
        ),
      tibble(AOD = rep(AOD_mid, length(AOD_hist))) %>%
        bind_cols(
          as_tibble(predict(
            AOD_ns3,
            rep(AOD_mid, length(AOD_hist))
          )) %>%
            set_names(~ str_c("AOD_ns3_", .x)),
          as_tibble(predict(
            AOD_ns5,
            rep(AOD_mid, length(AOD_hist))
          )) %>%
            set_names(~ str_c("AOD_ns5_", .x))
        )
    ) %>%
      mutate(AOD_x = AOD_hist)
  }, .id = "model") %>% 
    mutate(across(c(fit, lwr, upr), ~ expm1(.x) * 1e2))
  
  AOD_lab <- map_dfr(models, function(amodel) {
    relpred(
      amodel,
      tibble(AOD = rep(AOD_range[2], 10)) %>%
        bind_cols(
          as_tibble(predict(AOD_ns3, rep(AOD_range[2], 10))) %>%
            set_names(~ str_c("AOD_ns3_", .x)),
          as_tibble(predict(AOD_ns5, rep(AOD_range[2], 10))) %>%
            set_names(~ str_c("AOD_ns5_", .x))
        ),
      tibble(AOD = rep(AOD_mid, length(rep(AOD_range[2], 10)))) %>%
        bind_cols(
          as_tibble(predict(
            AOD_ns3,
            rep(AOD_mid, length(rep(AOD_range[2], 10)))
          )) %>%
            set_names(~ str_c("AOD_ns3_", .x)),
          as_tibble(predict(
            AOD_ns5,
            rep(AOD_mid, length(rep(AOD_range[2], 10)))
          )) %>%
            set_names(~ str_c("AOD_ns5_", .x))
        )
    ) %>%
      mutate(AOD_x = rep(AOD_range[2], 10))
  }, .id = "model") %>%
    group_by(model) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    separate(model, c("label", NA), "\\s", remove = F) %>% 
    mutate(across(c(fit, lwr, upr), ~ expm1(.x) * 1e2))
  
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
                        data = data_modeling,
                        bins = 60, fill = colors[2], color = "black", size = 0.3,
                        inherit.aes = F
    ) +
    ggside(x.pos = "bottom", collapse = "x") +
    scale_x_continuous(limits = AOD_range, expand = expansion(mult = c(0, 0.15))) +
    scale_xsidey_continuous(labels = NULL) +
    labs(subtitle = acrop)
  
  return(lst(p1, p2, plot_spacer()))
}

data <- data %>%
  mutate(plots = map2(fdata, crop_parent, get_plot))

patch <- reduce(flatten(data$plots), `+`) +
  plot_layout(ncol = 3, widths = c(1, 1, 0.7)) &
  scale_color_manual(
    values = wesanderson::wes_palette("Rushmore1",
                                      n = length(model_names),
                                      type = c("continuous")
    ),
    guide = guide_legend(reverse = T)
  ) &
  theme(plot.subtitle = element_text(size = 20))

ggsave("figures/response_function_Chinaplain.pdf", patch, width = 3, height = 2, scale = 5)
