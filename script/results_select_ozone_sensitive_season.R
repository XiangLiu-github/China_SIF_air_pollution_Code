source("script/loadPackages.R")
source("script/loadFunctions.R")
source("script/loadFormulas.R")

data <- get_data()

# must large than 3 months
data <- data %>%
  mutate(
    across(num_range("W126_", 4:7), ~ .x - W126_1, .names = "{.col}_m"),
    across(num_range("AOT40_", 4:7), ~ .x - AOT40_1, .names = "{.col}_m"),
    across(c(cloud, AOD), ~ .x^2, .names = "{.col}2")
  )

data %>%
  ggplot(aes(x = crop_parent, y = AOT40_7_m)) +
  geom_violin()

varss <-
  c(
    "GOSIF_sum", "cloud", "cloud2", "AOD", "AOD2",
    str_c("AOT40_", 3:7), str_c("AOT40_", 4:7, "_m"),
    str_c("bin", 0:39), str_c("surface_", 3:9)
  )

# because surface_1:2 is constant for crop 2
deFE <- data %>%
  mutate(surface_3 = surface_1 + surface_2 + surface_3) %>%
  drop_na(all_of(varss)) %>%
  nest(fdata = -c(crop, crop_parent)) %>%
  mutate(audata = map(fdata, ~ select(.x, fraction, x, y)))

deFE$fdata <- map2(deFE$fdata, deFE$audata, function(adata, aaudata) {
  pro_map_dfc(varss, function(aname) {
    fml <- str_c(aname, " ~ 1 | x_y[year]") %>% as.formula()

    tibble(!!aname := feols(fml, adata, weights = ~fraction)$residuals)
  }) %>%
    bind_cols(aaudata)
})

terms <- c(str_c("AOT40_", 3:7), str_c("AOT40_", 4:7, "_m"))
fmls <- map(terms, ~ str_c("GOSIF_sum ~ AOD + AOD ^ 2 + ", .x, " + cloud + cloud ^ 2 + ", add_terms(c(str_c("bin", 0:39), str_c("surface_", 3:9)))) %>%
  as.formula())

get_fdata <- function(adata) {
  adata <- adata %>% mutate(x_y = str_c(x, "_", y))

  set.seed(2022)
  folds <- group_vfold_cv(adata, group = x_y, v = 10)

  # ~ 2 mins
  results <-
    pro_map_dfr(folds$splits, function(asplit) {
      analy <- analysis(asplit)
      assess <- assessment(asplit)

      models <- map(fmls, ~ feols(.x, analy, lean = T, weights = ~fraction, nthreads = 0, notes = F)) %>%
        `names<-`(terms)

      # for each model
      map_dfr(models, function(amodel) {
        predicted <- predict(amodel, assess)

        return(predicted)
      }, .id = "model") %>%
        mutate(
          raw = assess$GOSIF_sum,
          fraction = assess$fraction
        )
    }, .id = "folds")

  gc()

  results %>%
    summarise(across(starts_with("AOT40"), ~ get_wr2_ss(raw, .x, fraction)))
}

deFE <- deFE %>%
  select(-audata) %>%
  unnest() %>%
  nest(fdata = -crop_parent) %>%
  mutate(data = map(fdata, get_fdata))

patch <- map(
  c("Maize", "Rice", "Wheat"),
  function(acrop) {
    deFE %>%
      filter(crop_parent == acrop) %>%
      select(-fdata) %>%
      unnest() %>%
      pivot_longer(-crop_parent, names_to = "model", values_to = "R2") %>%
      ggplot(aes(x = fct_reorder(model, R2), y = R2)) +
      geom_point() +
      theme_half_open(16, font_family = "Roboto Condensed") +
      background_grid() +
      xlab(NULL) +
      ylab(TeX("R$^{2}$")) +
      labs(subtitle = acrop)
  }
) %>%
  wrap_plots(ncol = 1)

ggsave("figures/ozone_sensitive_season.pdf", patch, width = 2, height = 1, scale = 7)









fml_base <- str_c("GOSIF_sum ~
            AOD + AOD ^ 2 +
            W126_7_m +
            cloud + cloud ^ 2 +
            ", add_terms(c(
  str_c("bin", 0:39),
  str_c("surface_", 1:9)
)), "| x_y[year]") %>%
  as.formula()

feols(fml_base, data, split = ~crop_parent, weights = ~fraction)


feols(GOSIF_sum ~ AOT40 | x_y[year], data, split = ~crop_parent, weights = ~fraction)[[1]]

a <- feols(GOSIF_sum ~ 1 | x_y[year], data, split = ~crop_parent, weights = ~fraction)[[1]]$residuals

b <- feols(AOT40 ~ 1 | x_y[year], data, split = ~crop_parent, weights = ~fraction)[[1]]$residuals

test <-
  feols(a ~ b, data = tibble(a = a, b = b, fraction = data %>% filter(crop_parent == "Maize") %>% pull(fraction)), weights = ~fraction)

rsq_vec(test$fitted.values, a)
get_wr2_ss(a, test$fitted.values, data %>% filter(crop_parent == "Maize") %>% pull(fraction))
