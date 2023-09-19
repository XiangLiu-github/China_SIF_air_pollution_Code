source("script/loadPackages.R")
source("script/loadFunctions.R")

data <- get_data()

data <- data %>%
  select(region, crop, x, y, year, `GR&EM`, MA) %>%
  mutate(gs = map2(`GR&EM`, MA, seq), .keep = "unused") %>%
  unnest()

p1 <- data %>%
  ggplot(aes(x = factor(gs), fill = region)) +
  facet_wrap(~crop) +
  geom_histogram(stat = "count", position = "stack", width = 0.7, color = "black", show.legend = F) +
  scale_x_discrete(name = "Month") +
  scale_y_continuous(name = "Count", labels = label_number()) +
  theme_half_open(18, font_family = "Roboto Condensed") +
  background_grid() +
  theme(legend.title = element_blank())

region <- read_xlsx("regions.xlsx") %>% mutate(region = str_c("R", region))

province_plot <- st_read("https://geo.datav.aliyun.com/areas_v3/bound/100000_full.json")

p2 <- inner_join(region, province_plot, by = c("province" = "name")) %>%
  ggplot(aes(geometry = geometry)) +
  facet_wrap(~crop_parent, nrow = 1) +
  geom_sf(data = province_plot) +
  geom_sf(aes(fill = region)) +
  theme_half_open(18, font_family = "Roboto Condensed") +
  background_grid(major = "none") +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 15),
    legend.position = c(1.1, 0.5),
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

patch <- p1 / p2 +
  plot_layout(guides = "collect") & 
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 40, family = "Roboto Condensed"))

ggsave("figures/growing_season.pdf", patch, width = 3, height = 3, scale = 5)
