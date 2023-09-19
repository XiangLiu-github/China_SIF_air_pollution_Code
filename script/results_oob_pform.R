source("script/loadPackages.R")
source("script/loadFunctions.R")

data <- get_data() %>%
  mutate(across(c(cloud, AOD, maxtmp, surface, root, prep),
    ~ .x^2,
    .names = "{.col}2"
  ))

varss <-
  c(
    "GOSIF_sum", "maxtmp", "maxtmp2", "surface", "surface2", "prep", "prep2",
    "root", "root2", "cloud", "cloud2", "AOD", "AOD2", "AOT40",
    str_c("bin", 0:39), str_c("step", 1:10),
    str_c("root_", 1:9), str_c("surface_", 3:9)
  )

# because surface_1:2 is constant for rice
deFE <- data %>%
  mutate(surface_3 = surface_1 + surface_2 + surface_3) %>%
  drop_na(all_of(varss)) %>%
  nest(fdata = -c(crop, crop_parent)) %>%
  mutate(audata = map(fdata, ~ select(.x, fraction, x, y)))

deFE$fdata <- pmap(list(deFE$crop_parent, deFE$fdata, deFE$audata), function(acrop, adata, aaudata) {
  map_dfc(varss, function(aname) {
    fml <- str_c(aname, " ~ 1 | x_y[year]") %>% as.formula()

    tibble(!!aname := feols(fml, adata, weights = ~fraction)$residuals)
  }, .progress = T) %>%
    bind_cols(aaudata)
})

fml_terms <- expand_grid(
  temperature = c(
    1, "maxtmp + maxtmp2",
    add_terms(str_c("bin", 0:39)), add_terms(str_c("step", 1:10))
  ),
  moisture = c(
    1, "prep + prep2", "surface + surface2", "root + root2",
    add_terms(str_c("root_", 1:9)), add_terms(str_c("surface_", 3:9))
  ),
  cloud = c(1, "cloud + cloud2"),
  aerosol = c(1, "AOD + AOD2"),
  ozone = c(1, "AOT40")
) %>%
  filter(!(temperature == "1" & moisture == "1" & cloud == "1" & aerosol == "1" & ozone == "1"))

fmls <- pmap_chr(fml_terms, ~ str_c("GOSIF_sum ~ ", add_terms(c(..1, ..2, ..3, ..4, ..5)))) %>%
  map(as.formula)

fml_trees <- fml_terms %>%
  array_tree()

get_datap <- function(adata) {
  adata <- adata %>% mutate(x_y = str_c(x, "_", y))

  set.seed(2022)
  folds <- group_vfold_cv(adata, group = x_y, v = 10)

  # ~ 8 mins
  results <-
    map_dfr(folds$splits, function(asplit) {
      analy <- analysis(asplit)
      assess <- assessment(asplit)

      models <- map(fmls, ~ feols(.x, analy,
        lean = T, weights = ~fraction,
        nthreads = 0, notes = F
      ))

      # for each model
      map2_dfr(models, fml_trees, function(amodel, atrees) {
        predicted <- predict(amodel, assess)

        tibble(
          R2 = get_wr2_ss(assess$GOSIF_sum, predicted, assess$fraction),
          name = list(atrees)
        )
      }, .id = "model")
    }, .id = "folds", .progress = T)

  gc()

  return(results)
}

deFE <- deFE %>%
  select(-audata) %>%
  unnest() %>%
  nest(fdata = -crop_parent) %>%
  mutate(plot_data = map(fdata, get_datap))

get_plot <- function(adata, acrop) {
  results_plot <- adata %>%
    unnest() %>%
    unnest() %>%
    filter(name != "1") %>%
    mutate(name = fcase(
      name == "maxtmp + maxtmp2", "temperature mean",
      name == add_terms(str_c("bin", 0:39)), "temperature 1 °C bin",
      name == add_terms(str_c("step", 1:10)), "temperature 4 °C bin",
      name == "prep + prep2", "prepcipitation",
      name == add_terms(str_c("root_", 1:9)), "root soil moisture bin",
      name == "root + root2", "root soil moisture mean",
      name == add_terms(str_c("surface_", 3:9)), "surface soil moisture bin",
      name == "surface + surface2", "surface soil moisture mean",
      name == "cloud + cloud2", "cloud",
      name == "AOD + AOD2", "aerosol",
      name == "AOT40", "ozone"
    )) %>%
    group_by(folds, model) %>%
    summarize(across(c(RMSE, R2), mean), name = add_terms(name), .groups = "drop")

  order_plot <- results_plot %>%
    group_by(model, name) %>%
    summarize(across(c(RMSE, R2), mean), .groups = "drop") %>%
    arrange(R2) %>%
    pull(name)

  results_plot <- results_plot %>%
    mutate(name = factor(name, levels = order_plot))

  results_plot %>%
    ggplot(aes(x = name, y = R2, fill = name)) +
    stat_summary(
      fun = "mean", geom = "bar", color = "black",
      show.legend = F, width = 0.7
    ) +
    # stat_summary(aes(label = sprintf("%0.4f", stat(y))),
    #              fun = function(x) round(mean(x), 4),
    #              geom = "text", family = "Roboto Condensed", size = 1,
    #              show.legend = F,
    #              position = position_nudge(y = 0.03)
    # ) +
    stat_summary(
      geom = "linerange", fun.max = "max", fun.min = "min",
      size = 1
    ) +
    axis_combmatrix(
      sep = " \\+ ",
      levels = c(
        "prepcipitation",
        "surface soil moisture bin", "surface soil moisture mean",
        "root soil moisture bin", "root soil moisture mean",
        "temperature mean", "temperature 1 °C bin", "temperature 4 °C bin",
        "cloud", "aerosol", "ozone"
      )
    ) +
    theme_half_open(25, font_family = "Roboto Condensed") +
    scale_y_continuous(expand = expansion()) +
    background_grid() +
    ylab(TeX("out-of-sample within R$^{2}$")) +
    xlab(NULL) +
    theme_combmatrix(
      combmatrix.label.text = element_text(
        family = "Roboto Condensed",
        size = 20
      ),
      combmatrix.panel.point.color.fill = "black",
      combmatrix.panel.line.size = 0,
      combmatrix.panel.point.size = 4
    ) +
    scale_fill_manual(
      values =
        NatParksPalettes::natparks.pals(
          "Arches",
          length(unique(results_plot$name))
        )
    ) +
    labs(tag = acrop)
}

deFE <- deFE %>%
  mutate(plot = map2(plot_data, crop_parent, get_plot))

p <- reduce(deFE$plot, `+`) +
  plot_layout(ncol = 1) &
  theme(plot.tag = element_text(size = 50))

ggsave("figures/oob_performance.pdf", p, width = 9, height = 5, scale = 5)


get_plot2 <- function(adata, acrop) {
  st_model <- "104"
  or_model <- c("128", "136")
  pc_model <- "62"
  xl_model <- "59"
  our_model <- "143"


  orderx <-
    inner_join(
      rownames_to_column(fml_terms, "model"),
      adata
    ) %>%
    group_by(model) %>%
    summarise(R2 = mean(R2)) %>%
    arrange(R2) %>%
    pull(model)

  ordery <- c(
    "prepcipitation",
    "surface soil moisture bin", "surface soil moisture mean",
    "root soil moisture bin", "root soil moisture mean",
    "temperature mean", "temperature 1 °C bin", "temperature 4 °C bin",
    "cloud", "aerosol", "ozone"
  )

  p1 <- inner_join(
    rownames_to_column(fml_terms, "model"),
    adata
  ) %>%
    mutate(model = factor(model, orderx)) %>%
    ggplot(aes(x = model, y = R2, fill = model)) +
    stat_summary(
      fun = "mean", geom = "bar", color = "black",
      show.legend = F, width = 0.7
    ) +
    stat_summary(
      geom = "linerange", fun.max = "max", fun.min = "min",
      size = 1
    ) +
    scale_fill_manual(
      values =
        NatParksPalettes::natparks.pals(
          "Arches",
          length(orderx)
        )
    ) +
    theme_half_open(25, font_family = "Roboto Condensed") +
    background_grid() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank()
    ) +
    scale_y_continuous(expand = expansion()) +
    ylab(TeX("out-of-sample within R$^{2}$")) +
    labs(tag = acrop)

  p2 <-
    fml_terms %>%
    rownames_to_column("model") %>%
    pivot_longer(-model) %>%
    mutate(model = factor(model, orderx)) %>%
    mutate(name = fcase(
      value == "maxtmp + maxtmp2", "temperature mean",
      value == add_terms(str_c("bin", 0:39)), "temperature 1 °C bin",
      value == add_terms(str_c("step", 1:10)), "temperature 4 °C bin",
      value == "prep + prep2", "prepcipitation",
      value == add_terms(str_c("root_", 1:9)), "root soil moisture bin",
      value == "root + root2", "root soil moisture mean",
      value == add_terms(str_c("surface_", 3:9)), "surface soil moisture bin",
      value == "surface + surface2", "surface soil moisture mean",
      value == "cloud + cloud2", "cloud",
      value == "AOD + AOD2", "aerosol",
      value == "AOT40", "ozone"
    )) %>%
    drop_na() %>%
    pivot_wider(values_fill = "1") %>%
    pivot_longer(-model) %>%
    mutate(
      color = fifelse(value == "1", "1", "0"),
      color = fifelse((model %in% st_model & color == "0"), "2", color),
      color = fifelse((model %in% or_model & color == "0"), "3", color),
      color = fifelse((model %in% pc_model & color == "0"), "4", color),
      color = fifelse((model %in% xl_model & color == "0"), "5", color),
      color = fifelse((model %in% our_model & color == "0"), "6", color),
      name = factor(name, levels = rev(ordery))
    ) %>%
    ggplot(aes(y = name, x = model, color = color)) +
    geom_point(size = 5, show.legend = F) +
    scale_color_manual(values = c("black", "grey", "green", "yellow", "purple", "blue", "red")) +
    theme_half_open(25, font_family = "Roboto Condensed") +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank()
    ) +
    geom_rect(
      aes(xmin = 0, xmax = length(orderx) + 0.5, ymin = ystart, ymax = yend, fill = col),
      data = data.frame(
        ystart = seq(0.5, length(ordery) - 0.5, 1),
        yend = seq(1.5, length(ordery) + 0.5, 1),
        col = rep_len(c("0", "1"), length(ordery))
      ),
      show.legend = F, inherit.aes = F,
      alpha = 0.1
    ) +
    scale_fill_manual(values = c("grey", "white"))

  return(lst(p1, p2))
}

deFE <- deFE %>%
  mutate(plot2 = map2(plot_data, crop_parent, get_plot2))

p <- reduce(flatten(deFE$plot2), `+`) +
  plot_layout(ncol = 1) &
  theme(plot.tag = element_text(size = 30))

ggsave("figures/oob_performance.pdf", p, width = 9, height = 5, scale = 5)
