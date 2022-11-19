source("script/loadPackages.R")
source("script/loadFunctions.R")
source("script/loadFormulas.R")

data <- get_data()

n <- 500

data <- data %>% nest(fadata = -crop_parent)

data$results <- pro_map(data$fadata, function(adata) {
  temp_year <- adata %>%
    select(year, GOSIF_sum) %>%
    nest_by(year) %>%
    ungroup()

  temp_x_y <- adata %>%
    select(x_y, GOSIF_sum) %>%
    nest_by(x_y) %>%
    ungroup()

  year_results <- map_dfr(1:n, function(aseed) {
    set.seed(aseed)
    rdm_year <- temp_year %>%
      sample_frac(1) %>%
      unnest()

    model_data <- adata %>%
      select(-c(year, GOSIF_sum)) %>%
      bind_cols(rdm_year)

    print(aseed)

    feols(fml_base, data = model_data, weights = ~fraction, lean = T, nthreads = 0, notes = F) %>%
      tidy()
  }, .id = "id")

  location_results <- map_dfr(1:n, function(aseed) {
    set.seed(aseed)
    rdm_x_y <- temp_x_y %>%
      sample_frac(1) %>%
      unnest()

    model_data <- adata %>%
      select(-c(x_y, GOSIF_sum)) %>%
      bind_cols(rdm_x_y)

    print(aseed)

    feols(fml_base, data = model_data, weights = ~fraction, lean = T, nthreads = 0, notes = F) %>%
      tidy()
  }, .id = "id")

  lst(year_results, location_results)
})

data <- data %>% select(-fadata)

saveRDS(data, "data/placebo_test.rds")

data <- read_rds("data/placebo_test.rds")

plot_hist <- function(adata) {
  adata %>%
    filter(str_detect(term, "AOD|AOT40")) %>%
    select(-c(std.error, statistic, p.value)) %>%
    pivot_wider(names_from = term, values_from = estimate) %>%
    ggplot(aes(x = AOT40 * 1e8)) +
    facet_wrap(~crop_parent, scales = "free") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_histogram(aes(fill = after_stat(..ncount..)), bins = 50, show.legend = F) +
    xlab(latex2exp::TeX("AOT40 coefficient ($\\times 10 ^{-8}$)")) +
    scale_fill_gradientn(colors = RColorBrewer::brewer.pal(7, "Blues")) +
    scale_y_continuous(expand = expansion()) +
    scale_x_continuous(limits = symmetric_limits)
}

plot_point <- function(adata) {
  adata %>%
    filter(str_detect(term, "AOD|PM25|AOT40")) %>%
    select(-c(std.error, statistic, p.value)) %>%
    pivot_wider(names_from = term, values_from = estimate) %>%
    ggplot(aes(x = AOD * 1e3, y = `I(AOD^2)` * 1e3)) +
    facet_wrap(~crop_parent, scales = "free") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_hdr(aes(fill = after_stat(probs)), show.legend = F) +
    geom_hdr_points(aes(color = after_stat(probs)), show.legend = F) +
    geom_hdr_rug(aes(fill = after_stat(probs)), show.legend = F) +
    xlab(latex2exp::TeX("AOD coefficient ($\\times 10 ^{-3}$)")) +
    ylab(latex2exp::TeX("AOD$^{2}$ coefficient ($\\times 10 ^{-3}$)")) +
    scale_fill_brewer(palette = "Greens") +
    scale_color_brewer(palette = "Greens") +
    scale_x_continuous(limits = symmetric_limits) +
    scale_y_continuous(limits = symmetric_limits)
}

p1 <-
  data %>%
  unnest_wider(results) %>%
  select(crop_parent, year_results) %>%
  unnest() %>%
  plot_hist() +
  labs(subtitle = "reshuffled by years")

p2 <-
  data %>%
  unnest_wider(results) %>%
  select(crop_parent, year_results) %>%
  unnest() %>%
  plot_point() +
  labs(subtitle = "reshuffled by years")

p3 <-
  data %>%
  unnest_wider(results) %>%
  select(crop_parent, location_results) %>%
  unnest() %>%
  plot_hist() +
  labs(subtitle = "reshuffled by locations")

p4 <-
  data %>%
  unnest_wider(results) %>%
  select(crop_parent, location_results) %>%
  unnest() %>%
  plot_point() +
  labs(subtitle = "reshuffled by locations")

patch <- p1 / p3 / p2 / p4 &
  theme_half_open(18, font_family = "Roboto Condensed") &
  background_grid() &
  theme(plot.subtitle = element_text(size = 20))

ggsave("figures/placebo_test.pdf", patch, width = 3, height = 4, scale = 4)
