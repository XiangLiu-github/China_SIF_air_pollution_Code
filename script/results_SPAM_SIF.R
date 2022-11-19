source("script/loadPackages.R")
source("script/loadFunctions.R")

data <- qread("data/tidied.qs", nthreads = qn) %>%
  mutate(crop = fifelse(str_detect(crop, "Rice"), "Rice", crop)) %>%
  rename_with(~ str_c(.x, "_mean"), c(GOSIF, RTSIF, Wenetal, CSIF))

SIF <- data %>%
  lazy_dt() %>%
  group_by(x, y, year, crop, county, county_code, province_code) %>%
  summarise(across(c(starts_with("GOSIF"), fraction), mean), .groups = "drop") %>%
  as_tibble()

spam_all <-
  map2_dfr(
    c("Wheat", "Rice", "Maize"),
    c("WHEA", "RICE", "MAIZ"),
    function(acrop1, acrop) {
      mask <- rast("data/outputs/masks/mask_extraction.tif")[[acrop1]]

      harv <- c(
        str_c("/vsizip/../data_archive/SPAM/spam2000v3.0.7_global_harvested-area.geotiff.zip/spam2000V3r107_global_H_", acrop, "_A.tif"),
        str_c("/vsizip/../data_archive/SPAM/spam2005v3r2_global_harv_area.geotiff.zip/geotiff_global_harv_area/SPAM2005V3r2_global_H_TA_", acrop, "_A.tif"),
        str_c("/vsizip/../data_archive/SPAM/spam2010v2r0_global_harv_area.geotiff.zip/spam2010V2r0_global_H_", acrop, "_A.tif")
      ) %>%
        map(rast) %>%
        map(~ project(.x, mask)) %>%
        rast() %>%
        mask(mask) %>%
        `names<-`(seq(2000, 2010, by = 5)) %>%
        as.data.frame(xy = T, na.rm = F) %>%
        pivot_longer(-c(x, y),
          names_to = "year", names_transform = list(year = as.integer),
          values_to = "area", values_drop_na = T
        ) %>%
        filter(area != 0)

      prod <- c(
        str_c("/vsizip/../data_archive/SPAM/spam2000v3.0.7_global_production.geotiff.zip/spam2000V3r107_global_R_", acrop, "_A.tif"),
        str_c("/vsizip/../data_archive/SPAM/spam2005v3r2_global_prod.geotiff.zip/geotiff_global_prod/SPAM2005V3r2_global_P_TA_", acrop, "_A.tif"),
        str_c("/vsizip/../data_archive/SPAM/spam2010v2r0_global_prod.geotiff.zip/spam2010V2r0_global_P_", acrop, "_A.tif")
      ) %>%
        map(rast) %>%
        map(~ project(.x, mask)) %>%
        rast() %>%
        mask(mask) %>%
        `names<-`(seq(2000, 2010, by = 5)) %>%
        as.data.frame(xy = T, na.rm = F) %>%
        pivot_longer(-c(x, y),
          names_to = "year", names_transform = list(year = as.integer),
          values_to = "production", values_drop_na = T
        ) %>%
        filter(production != 0)

      inner_join(harv, prod) %>%
        mutate(yield = production / area, crop = acrop1)
    }
  ) %>%
  trim_xy()

joined <- inner_join(SIF, spam_all) %>%
  drop_na() %>%
  group_by(crop) %>%
  filter(yield < mean(yield) + 3 * sd(yield) & yield > mean(yield) - 3 * sd(yield)) %>%
  ungroup()

hist(joined$yield[joined$crop == "Maize"])

p1 <-
  filter(joined, fraction >= max(joined$fraction) * 0.7) %>%
  pivot_longer(contains("GOSIF"), names_transform = list(name = ~ str_replace(.x, "_", " "))) %>%
  ggplot(aes(y = yield, x = value)) +
  facet_grid2(vars(name), vars(crop), scales = "free", independent = "all") +
  geom_pointdensity(aes(color = stat(ndensity))) +
  geom_smooth(method = "lm", se = F) +
  stat_poly_eq(use_label(c("n", "R2", "p.value")), show.legend = F, family = "Roboto Condensed") +
  theme_half_open(18, font_family = "Roboto Condensed") +
  background_grid() +
  scale_color_gradientn(colours = PNWColors::pnw_palette("Moth") %>% rev()) +
  xlab(TeX("SIF (mW $m^{-2}$ $nm^{-1}$ $sr^{-1}$)")) +
  ylab(TeX("Yield (kg/ha)"))

ggsave("figures/spam_validataion.pdf", p1, width = 1, height = 1, scale = 10)
