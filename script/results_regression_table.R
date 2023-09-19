source("script/loadPackages.R")
source("script/loadFunctions.R")
source("script/loadFormulas.R")

data <- get_data()

model <-
  feols(fml_base,
    data,
    cluster = ~county,
    split = ~crop_parent,
    weights = ~fraction,
    nthreads = 0,
    lean = T
  )

etable(model,
  keep = "AOD|AOT40", tex = F,
  dict = c(
    "GOSIF_sum" = "ln(GOSIF)",
    "crop_parent" = "crop"
  )
)

etable(model,
  keep = "AOD|AOT40", tex = TRUE,
  dict = c(
    "GOSIF_sum" = "ln(GOSIF)",
    "crop_parent" = "crop"
  )
)

model <-
  feols(fml_base_w126,
    data,
    cluster = ~county,
    split = ~crop_parent,
    weights = ~fraction,
    nthreads = 0,
    lean = T
  )

etable(model,
  keep = "AOD|W126", tex = TRUE,
  dict = c(
    "GOSIF_sum" = "ln(GOSIF)",
    "crop_parent" = "crop"
  )
)


model <-
  feols(fml_base,
    data %>% mutate(AOT40 = AOT40_7),
    cluster = ~county,
    split = ~crop_parent,
    weights = ~fraction,
    nthreads = 0,
    lean = T
  )

etable(model,
  keep = "AOD|AOT40", tex = F,
  dict = c(
    "GOSIF_sum" = "ln(GOSIF)",
    "crop_parent" = "crop"
  )
)

etable(model,
  keep = "AOD|AOT40", tex = TRUE,
  dict = c(
    "GOSIF_sum" = "ln(GOSIF)",
    "crop_parent" = "crop"
  )
)

# data_irg <- data %>%
#   filter(crop_parent != "Rice") %>%
#   mutate(irg = fcase(irg_fraction < 0.1, "rainfed", irg_fraction > 0.75, "irrigated")) %>%
#   drop_na(irg) %>%
#   mutate(class = str_c(crop_parent, "_", irg))
#
# data_irg %>%
#   ggplot(aes(x = irg, y = fraction)) +
#   facet_wrap(~crop_parent) +
#   geom_violin()
#
# model_irg <-
#   feols(fml_base,
#     data_irg,
#     cluster = ~county,
#     split = ~class,
#     weights = ~fraction,
#     nthreads = 0,
#     lean = T
#   )
#
# etable(model_irg,
#   keep = "AOD|AOT40", tex = TRUE,
#   dict = c(
#     "GOSIF_sum" = "ln(GOSIF)",
#     "crop_parent" = "crop"
#   )
# )
#
#
# data_sm <- data %>%
#   group_by(x_y) %>%
#   mutate(surface_mean = mean(surface)) %>%
#   group_by(crop_parent) %>%
#   mutate(surface_mean = fcase(
#     surface_mean >= quantile(surface, 0.75), "wet",
#     surface_mean <= quantile(surface, 0.25), "dry"
#   )) %>%
#   ungroup() %>%
#   drop_na(surface_mean) %>%
#   mutate(class = str_c(crop_parent, "_", surface_mean))
#
# data_sm %>%
#   ggplot(aes(x = surface_mean, y = fraction)) +
#   facet_wrap(~crop_parent) +
#   geom_violin()
#
# model_sm <-
#   feols(fml_base,
#     data_sm,
#     cluster = ~county,
#     split = ~class,
#     weights = ~fraction,
#     nthreads = 0,
#     lean = T
#   )
#
# etable(model_sm,
#   keep = "AOD|AOT40", tex = TRUE,
#   dict = c(
#     "GOSIF_sum" = "ln(GOSIF)",
#     "crop_parent" = "crop"
#   )
# )
#
#
# data_temp <- data %>%
#   group_by(x_y) %>%
#   mutate(maxtmp_mean = mean(maxtmp)) %>%
#   group_by(crop_parent) %>%
#   mutate(maxtmp_mean = fcase(
#     maxtmp_mean >= quantile(maxtmp, 0.75), "hot",
#     maxtmp_mean <= quantile(maxtmp, 0.25), "cold"
#   )) %>%
#   ungroup() %>%
#   drop_na(maxtmp_mean) %>%
#   mutate(class = str_c(crop_parent, "_", maxtmp_mean))
#
# data_temp %>%
#   ggplot(aes(x = maxtmp_mean, y = fraction)) +
#   facet_wrap(~crop_parent) +
#   geom_violin()
#
# model_temp <-
#   feols(fml_base,
#     data_temp,
#     cluster = ~county,
#     split = ~class,
#     weights = ~fraction,
#     nthreads = 0,
#     lean = T
#   )
#
# etable(model_temp,
#   keep = "AOD|AOT40", tex = TRUE,
#   dict = c(
#     "GOSIF_sum" = "ln(GOSIF)",
#     "crop_parent" = "crop"
#   )
# )
