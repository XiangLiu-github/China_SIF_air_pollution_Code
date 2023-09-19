source("script/loadPackages.R")
source("script/loadFunctions.R")

shp <- st_read("data/inputs/shp/2020省矢量.shp")

site_info <- read_csv("../data_archive/GEBA/GEBA_2023-03-28_10-01-07_metadata.csv") %>%
  select(tskey, ebcode, sgxlon, sgylat) %>%
  mutate(ebcode = case_when(
    ebcode == "GLOBAL" ~ "global",
    ebcode == "DIFFUS" ~ "diffuse",
    TRUE ~ NA
  )) %>%
  drop_na()

sites <- st_as_sf(site_info, coords = c("sgxlon", "sgylat"), remove = F)

site_data <- read_csv("../data_archive/GEBA/GEBA_2023-03-28_10-01-07_monthly_data.csv") %>%
  filter(computed_flag_avg == 8) %>%
  select(-computed_flag_avg) %>%
  rename(insitu = converted_flux_avg)

SR_data <-
  list.files("../data_archive/China_5km_radition/", full.names = T) %>%
  map(function(afile) {
    date <- str_remove(afile, "../data_archive/China_5km_radition/RAD_") %>%
      str_remove(".h5")

    data_mx <- rast(afile)[]
    data_mx <- data_mx * 0.01

    data_rst <- brick(
      nrows = 901, ncols = 1401,
      xmn = 71.025 - 0.05 / 2, xmx = 141.025 + 0.05 / 2,
      ymn = 14.975 - 0.05 / 2, ymx = 59.975 + 0.05 / 2,
      nl = 2
    )

    values(data_rst) <- data_mx
    names(data_rst) <- str_c(names(data_rst), "_", date)

    return(rast(data_rst))
  }, .progress = T)

SR_data <- rast(SR_data)

SR_data <- terra::extract(SR_data, sites, ID = F) %>%
  bind_cols(sites %>% st_drop_geometry()) %>%
  pivot_longer(starts_with("monthly"),
    names_to = c(NA, NA, "type", NA, "date"), names_sep = "_", names_transform = list(date = ~ ymd(.x, truncated = 4)),
    values_to = "old"
  ) %>%
  mutate(year = year(date), month = month(date), .keep = "unused") %>%
  filter(ebcode == type)

BaRAD <- list.files("../data_archive/BaRAD", full.names = T) %>%
  rast() %>%
  `names<-`(str_c(names(.), "_", rep(1980:2019, each = 3 * 12))) %>%
  terra::extract(sites, ID = F) %>%
  bind_cols(sites %>% st_drop_geometry()) %>%
  pivot_longer(starts_with("K_"),
    names_to = c(NA, "type", "month", "year"), names_sep = "_", names_transform = list(month = as.integer, year = as.integer),
    values_to = "new"
  ) %>%
  mutate(type = case_when(
    type == "down" ~ "global",
    type == "diff" ~ "diffuse",
    TRUE ~ NA
  )) %>%
  drop_na() %>%
  filter(ebcode == type)

inner_join(SR_data, site_data) %>%
  ggplot(aes(x = old, y = insitu)) +
  facet_grid2(vars(type), vars(month), scales = "free", independent = "all") +
  geom_abline(slope = 1, intercept = 0) +
  geom_point() +
  stat_poly_eq()

inner_join(BaRAD, site_data) %>%
  filter(year %in% 2007:2018) %>%
  ggplot(aes(x = new, y = insitu)) +
  facet_grid2(vars(type), vars(month), scales = "free", independent = "all") +
  geom_abline(slope = 1, intercept = 0) +
  geom_point() +
  stat_poly_eq()

a <- inner_join(SR_data, site_data) %>%
  filter(year %in% 2007:2018) %>%
  group_by(type, month) %>%
  rsq(old, insitu) %>%
  select(-c(.metric, .estimator)) %>%
  pivot_wider(names_from = type, names_prefix = "Jiang et al. ", values_from = .estimate)

b <- inner_join(BaRAD, site_data) %>%
  group_by(type, month) %>%
  rsq(new, insitu) %>%
  select(-c(.metric, .estimator)) %>%
  pivot_wider(names_from = type, names_prefix = "Chakraborty et al. ", values_from = .estimate)

inner_join(a, b) %>%
  kableExtra::kbl(
    format = "latex",
    digits = c(0, 3, 3, 3, 3)
  )
