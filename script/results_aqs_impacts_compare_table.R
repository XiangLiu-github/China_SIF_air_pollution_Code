source("script/loadPackages.R")
source("script/loadFunctions.R")

inner_join(
  read_rds("data/impacts_aerosol_summarised.rds") %>%
    unnest_wider(AOD_results) %>%
    select(-province_level) %>%
    unnest() %>%
    group_by(crop_parent, region, pm_level) %>%
    summarise(across(contains("%"), fmean), .groups = "drop") %>%
    filter(pm_level == 35) %>%
    select(crop_parent, region, `50%`) %>%
    rename(nonhete = `50%`),
  read_rds("data/impacts_aerosol_hete_full_summarised.rds") %>%
    unnest_wider(AOD_results) %>%
    select(-province_level) %>%
    unnest() %>%
    group_by(crop_parent, region, pm_level) %>%
    summarise(across(contains("%"), fmean), .groups = "drop") %>%
    filter(pm_level == 35) %>%
    select(crop_parent, region, `50%`) %>%
    rename(hete = `50%`)
) %>%
  mutate(across(c(nonhete, hete), ~ .x * 1e2)) %>%
  kableExtra::kbl(
    col.names = c("Crop", "Region", "Non-interactive model", "Interactive model"),
    format = "latex",
    digits = c(0, 0, 2, 2),
    toprule = "\\toprule",
    bottomrule = "\\bottomrule",
    midrule = "\\midrule"
  )

inner_join(
  read_rds("data/impacts_ozone_summarised.rds") %>%
    unnest_wider(ozone_results) %>%
    select(-province_level) %>%
    unnest() %>%
    group_by(crop_parent, region, peak_level) %>%
    summarise(across(contains("%"), fmean), .groups = "drop") %>%
    filter(peak_level == 60) %>%
    select(crop_parent, region, `50%`) %>%
    rename(nonhete = `50%`),
  read_rds("data/impacts_ozone_hete_full_summarised.rds") %>%
    unnest_wider(ozone_results) %>%
    select(-province_level) %>%
    unnest() %>%
    group_by(crop_parent, region, peak_level) %>%
    summarise(across(contains("%"), fmean), .groups = "drop") %>%
    filter(peak_level == 60) %>%
    select(crop_parent, region, `50%`) %>%
    rename(hete = `50%`)
) %>%
  mutate(across(c(nonhete, hete), ~ .x * 1e2)) %>%
  kableExtra::kbl(
    col.names = c("Crop", "Region", "Non-interactive model", "Interactive model"),
    format = "latex",
    digits = c(0, 0, 2, 2),
    toprule = "\\toprule",
    bottomrule = "\\bottomrule",
    midrule = "\\midrule"
  )
