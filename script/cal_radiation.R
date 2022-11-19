source("script/loadPackages.R")
source("script/loadFunctions.R")

mask <- rast("data/outputs/masks/mask_extraction.tif") %>% sum(na.rm = T)

SR_data <-
  list.files("D:/data/China_5km_radition/", full.names = T) %>%
  pro_map(function(afile) {
    date <- str_remove(afile, "D:/data/China_5km_radition/RAD_") %>%
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
  })

SR <- rast(SR_data)

DR <- SR["diffuse_radiation"] %>%
  resample(mask, method = "near") %>%
  mask(mask)
GR <- SR["global_radiation"] %>%
  resample(mask, method = "near") %>%
  mask(mask)

rm(SR_data, SR)

gc()

DR <- DR %>%
  as.data.frame(xy = T) %>%
  pivot_longer(-c(x, y),
    names_prefix = "monthly_total_diffuse_radiation_",
    names_to = "date", names_transform = list(date = ym),
    values_to = "DR"
  ) %>%
  mutate(year = year(date), month = month(date), .keep = "unused")

GR <- GR %>%
  as.data.frame(xy = T) %>%
  pivot_longer(-c(x, y),
    names_prefix = "monthly_total_global_radiation_",
    names_to = "date", names_transform = list(date = ym),
    values_to = "GR"
  ) %>%
  mutate(year = year(date), month = month(date), .keep = "unused")

SR <- inner_join(DR, GR) %>%
  trim_xy()

check_join(SR)

saveRDS(SR, "data/outputs/radiation/tidied.rds")
