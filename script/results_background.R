source("script/loadPackages.R")
source("script/loadFunctions.R")

order_crop <- c("Maize", "Rice(SR&ER)", "Rice(LR)", "Wheat")
SIF_color <- "BlueYl"
W126_color <- "GrPink2"
AOD_color <- "BlueOr"

province <- read_rds("data/inputs/shp/province.rds")
jiuduanxian <- read_rds("data/inputs/shp/jiuduanxian.rds")

data <- get_data()

test_data <- data %>%
  group_by(x, y, crop) %>%
  summarise(across(c(GOSIF, fraction, AOD, W126), ~ mean(.x, na.rm = T)),
    .groups = "drop"
  )

data_SIF <- test_data %>%
  nest(fdata = -crop) %>%
  mutate(fdata = map(fdata, ~ bi_class(.x,
    x = GOSIF, y = fraction,
    style = "fisher", dim = 4
  ))) %>%
  unnest() %>%
  mutate(crop = factor(crop, order_crop))

data_W126 <- test_data %>%
  nest(fdata = -crop) %>%
  mutate(fdata = map(fdata, ~ bi_class(.x,
    x = W126, y = fraction,
    style = "fisher", dim = 4
  ))) %>%
  unnest() %>%
  mutate(crop = factor(crop, order_crop))

data_AOD <- test_data %>%
  nest(fdata = -crop) %>%
  mutate(fdata = map(fdata, ~ bi_class(.x,
    x = AOD, y = fraction,
    style = "fisher", dim = 4
  ))) %>%
  unnest() %>%
  mutate(crop = factor(crop, order_crop))

map_theme <- function(p) {
  p +
    xlab(NULL) +
    ylab(NULL) +
    theme_ipsum_rc(
      panel_spacing = grid::unit(0, "lines"),
      plot_margin = margin(0, 0, 0, 0),
      grid = F, axis = F,
      base_size = 15,
      axis_title_size = 15,
      strip_text_size = 15
    ) +
    theme(
      axis.text.x = element_blank(), axis.text.y = element_blank()
    )
}

legend_theme <- function(p) {
  p +
    theme(
      axis.title = element_text(family = "Roboto Condensed", size = 5),
      axis.line.y = element_line(arrow = arrow(
        length = unit(0.1, "cm"),
        ends = "last", type = "open"
      )),
      axis.line.x = element_line(arrow = arrow(
        length = unit(0.1, "cm"),
        ends = "last", type = "open"
      ))
    )
}

map_SIF <- data_SIF %>%
  ggplot(aes(x = x, y = y, fill = bi_class)) +
  facet_wrap(~crop, nrow = 1) +
  geom_sf(data = province, inherit.aes = F, fill = NA) +
  geom_sf(data = jiuduanxian, inherit.aes = F, fill = NA) +
  geom_tile(show.legend = FALSE) +
  bi_scale_fill(pal = SIF_color, dim = 4)


legend_SIF <- bi_legend(
  pal = SIF_color,
  dim = 4,
  xlab = "Higher SIF",
  ylab = "Higher crop fraction",
  size = 8,
  arrows = F,
  pad_width = 0.5
)

map_W126 <- data_SIF %>%
  ggplot(aes(x = x, y = y, fill = bi_class)) +
  facet_wrap(~crop, nrow = 1) +
  geom_sf(data = province, inherit.aes = F, fill = NA) +
  geom_sf(data = jiuduanxian, inherit.aes = F, fill = NA) +
  geom_tile(show.legend = FALSE) +
  bi_scale_fill(pal = W126_color, dim = 4)


legend_W126 <- bi_legend(
  pal = W126_color,
  dim = 4,
  xlab = "Higher W126",
  ylab = "Higher crop fraction",
  size = 8,
  arrows = F,
  pad_width = 0.5
)

map_AOD <- data_AOD %>%
  ggplot(aes(x = x, y = y, fill = bi_class)) +
  facet_wrap(~crop, nrow = 1) +
  geom_sf(data = province, inherit.aes = F, fill = NA) +
  geom_sf(data = jiuduanxian, inherit.aes = F, fill = NA) +
  geom_tile(show.legend = FALSE) +
  bi_scale_fill(pal = AOD_color, dim = 4)


legend_AOD <- bi_legend(
  pal = AOD_color,
  dim = 4,
  xlab = "Higher AOD",
  ylab = "Higher crop fraction",
  size = 8,
  arrows = F,
  pad_width = 0.5
)


patch <-
  (map_theme(map_SIF) +
    inset_element(legend_theme(legend_SIF),
      left = 0, bottom = 0, right = 0.08, top = 0.5
    )) /
    (map_theme(map_W126) +
      inset_element(legend_theme(legend_W126),
        left = 0, bottom = 0, right = 0.08, top = 0.5
      )) /
    (map_theme(map_AOD) +
      inset_element(legend_theme(legend_AOD),
        left = 0, bottom = 0, right = 0.08, top = 0.5
      ))

ggsave("figures/background.pdf", patch, width = 4, height = 3, scale = 4)
