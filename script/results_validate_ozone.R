source("script/loadPackages.R")
source("script/loadFunctions.R")

ozone <-
  inner_join(
    read_rds("../Vegetation_ozone/data/inputs/station_raw/ozone/China.rds"),
    read_rds("../pre_NS2021/data/inputs/station_raw/ozone/China.rds")
  ) %>%
  mutate(O3 = rollingo3 * 1e3 * 2)

ozone_xy <- ozone %>%
  distinct(Longitude, Latitude) %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = "EPSG:4326", remove = F)

Zhou <- list.files("../data_archive/China_10km_ozone/Zhouetal2022/", full.names = T)[11:16] %>%
  map(~ rast(.x) %>% `names<-`(time(.))) %>%
  rast()

Zhou <- terra::extract(Zhou, ozone_xy, ID = F) %>%
  bind_cols(ozone_xy %>% st_drop_geometry()) %>%
  pivot_longer(-c(Latitude, Longitude),
    names_to = "date",
    names_transform = list(date = ymd),
    values_to = "Zhou"
  ) %>%
  unique() %>%
  mutate(year = year(date), month = month(date), .keep = "unused")

p <- inner_join(ozone, Zhou) %>%
  ggplot(aes(y = O3, x = Zhou)) +
  geom_abline(slope = 1, intercept = 0, color = "blue") +
  geom_pointdensity(aes(color = stat(ndensity))) +
  geom_smooth(method = "lm", color = "red") +
  scale_y_continuous(name = "Observed monthly MDA8") +
  scale_x_continuous(name = "Predicted monthly MDA8 in Zhou et al.") +
  stat_poly_eq(use_label(c("eq.label", "rr.label", "p.value.label")), family = "Roboto Condensed", size = 5) +
  theme_half_open(16, font_family = "Roboto Condensed") +
  background_grid()

ggsave("figures/ozone_validation.pdf", p, width = 1.2, height = 1, scale = 5)
