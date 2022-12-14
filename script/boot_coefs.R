source("script/loadPackages.R")
source("script/loadFunctions.R")
source("script/loadFormulas.R")

# f1 ----------------------------------------------------------------------

data <- get_data()
n <- 500

data <- data %>% tidytable::nest_by.(crop_parent, .key = "fdata")

data <- data %>%
  mutate(coefs = map(fdata, function(adata) {
    adata <- adata %>%
      tidytable::nest_by.(county, .key = "ffdata")

    pro_map_dfr(1:n, function(anum) {
      set.seed(anum)
      aadata <- adata %>%
        tidytable::slice_sample.(prop = 1, replace = T) %>%
        tidytable::unnest.(ffdata)

      feols(fml_base, aadata,
        weights = ~fraction,
        nthreads = 0, notes = F, lean = T
      ) %>%
        tidy()
    }, .id = "id")
  })) %>%
  as_tibble()

saveRDS(data %>% select(-fdata), "data/boots_f1.rds")

# f2 ----------------------------------------------------------------------

# data = read_rds('data/stats.rds')
#
# data = data %>%
#   mutate(crop_parent = fcase(crop == '玉米', 'Maize',
#                              crop == '稻谷', 'Rice',
#                              crop == '小麦', 'Wheat')) %>%
#   nest(fdata = - crop_parent)
#
# data = data %>%
#   mutate(coefs = pro_map(fdata, function(adata) {
#
#   adata = adata %>% bootstraps(500)
#
#   map_dfr(adata$splits, function(asplit) {
#
#     feols(production ~ GOSIF_mean | province_name + year,
#           data = analysis(asplit), weights = ~ area, notes = F) %>%
#       tidy()
#
#   }, .id = 'id')
#
# }))
#
# saveRDS(data %>% select(-fdata), 'data/boots_f2.rds')
