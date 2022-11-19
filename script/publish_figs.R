rm(list = ls())
library(tidyverse)

# remove files ------------------------------------------------------------

list.files("figures_paper/", full.names = T) %>%
  walk(file.remove)

list.files("figures_png/", full.names = T) %>%
  walk(file.remove)

list.files("figures_paper_png/", full.names = T) %>%
  walk(file.remove)

# main --------------------------------------------------------------------

c(
  "figures/response_function.pdf",
  "figures/hist_impact_map.pdf",
  "figures/hist_impact_fig3.pdf",
  "figures/hist_impact_fig4.pdf"
) %>%
  iwalk(function(afig, anum) {
    file.copy(afig, str_c("figures_paper/Fig", anum, ".pdf"), overwrite = T)
  })

# supplementary ----------------------------------------------------------------

c(
  "figures/NF2022_workflow.pdf",
  "figures/background_line_Maize.pdf",
  "figures/background_line_Rice.pdf",
  "figures/background_line_Wheat.pdf",
  "figures/spam_validataion.pdf",
  "figures/EC_validation.pdf",
  "figures/background_map_Maize.pdf",
  "figures/background_map_Rice.pdf",
  "figures/background_map_Wheat.pdf",
  "figures/gam_performance.pdf",
  "figures/oob_performance.pdf",
  "figures/placebo_test.pdf",
  "figures/response_radiation.pdf",
  "figures/Step_2_counterfactual_air_pollution_level.pdf",
  "figures/compare_sensitivity.pdf",
  "figures/response_hete_irrigation.pdf",
  "figures/response_hete_temp.pdf",
  "figures/response_hete_surface.pdf",
  "figures/hist_impact_fig4_Feng2022.pdf",
  "figures/response_function_Chinaplain.pdf",
  "figures/hist_impact_matrix_new.pdf",
  "figures/ozone_sensitive_season.pdf",
  "figures/response_dSIF.pdf",
  "figures/response_climate.pdf",
  "figures/response_FE.pdf",
  "figures/hist_impact_2013.pdf"
) %>%
  iwalk(function(afig, anum) {
    file.copy(afig, str_c("figures_paper/SIFig", anum, ".pdf"), overwrite = T)
  })

# check -------------------------------------------------------------------

stopifnot(length(list.files("figures_paper/")) == length(list.files("figures/")))


# convert to png ----------------------------------------------------------

library(magick)
list.files("figures/", full.names = F) %>%
  walk(function(aname) {
    image_read_pdf(str_c("figures/", aname), density = 300) %>%
      image_write(path = str_c("figures_png/", str_replace(aname, ".pdf", ".png")), format = "png")
  })

list.files("figures_paper/", full.names = F) %>%
  walk(function(aname) {
    image_read_pdf(str_c("figures_paper/", aname), density = 300) %>%
      image_write(path = str_c("figures_paper_png/", str_replace(aname, ".pdf", ".png")), format = "png")
  })
