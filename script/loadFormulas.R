# base --------------------------------------------------------------------

fml_base <- str_c("GOSIF_sum ~
            AOD + AOD ^ 2 + 
            AOT40 +
            cloud + cloud ^ 2 +
            ", add_terms(c(
  str_c("bin", 0:39),
  str_c("surface_", 1:9)
)), "|
            x_y[year]") %>%
  as.formula()

fml_base_rad <- str_c("GOSIF_sum ~
            GR + GR ^ 2 +
            DR + DR ^ 2 +
            AOT40 +
            ", add_terms(c(
  str_c("bin", 0:39),
  str_c("surface_", 1:9)
)), "|
            x_y[year]") %>%
  as.formula()

# function ----------------------------------------------------------------

fml_linear <- str_c("GOSIF_sum ~
            AOD +
            AOT40 +
            cloud + cloud ^2 +
             ", add_terms(c(
  str_c("bin", 0:39),
  str_c("surface_", 1:9)
)), "| x_y[year]") %>%
  as.formula()

fml_poly3 <- str_c("GOSIF_sum ~
            AOD + AOD ^ 2 + AOD ^ 3 +
            AOT40 + AOT40 ^2 + AOT40 ^3 +
            cloud + cloud ^2 +
             ", add_terms(c(
  str_c("bin", 0:39),
  str_c("surface_", 1:9)
)), "| x_y[year]") %>%
  as.formula()

fml_ns3 <- str_c("GOSIF_sum ~
            ", add_terms(str_c("AOD_ns3_", 1:3)), " +
            ", add_terms(str_c("AOT40_ns3_", 1:3)), " +
            cloud + cloud ^2 +
             ", add_terms(c(
  str_c("bin", 0:39),
  str_c("surface_", 1:9)
)), "| x_y[year]") %>%
  as.formula()

fml_ns5 <- str_c("GOSIF_sum ~
            ", add_terms(str_c("AOD_ns5_", 1:5)), " +
            ", add_terms(str_c("AOT40_ns5_", 1:5)), " +
            cloud + cloud ^2 +
             ", add_terms(c(
  str_c("bin", 0:39),
  str_c("surface_", 1:9)
)), "| x_y[year]") %>%
  as.formula()

fml_noaerosol <- str_c("GOSIF_sum ~
            AOT40 +
            cloud + cloud ^2 +
            ", add_terms(c(
  str_c("bin", 0:39),
  str_c("surface_", 1:9)
)), "| x_y[year]") %>%
  as.formula()

fml_noozone <- str_c("GOSIF_sum ~
            AOD + AOD ^ 2 +
            cloud + cloud ^2 +
            ", add_terms(c(
  str_c("bin", 0:39),
  str_c("surface_", 1:9)
)), "| x_y[year]") %>%
  as.formula()

# FE ----------------------------------------------------------------------

fml_grid_year_trend2 <- str_c("GOSIF_sum ~
            AOD + AOD ^ 2 +
            AOT40 +
            cloud + cloud ^ 2 +
            ", add_terms(c(
  str_c("bin", 0:39),
  str_c("surface_", 1:9)
)), "|
            x_y[year + year ^ 2]") %>%
  as.formula()

fml_no_year <- str_c("GOSIF_sum ~
            AOD + AOD ^ 2 +
            AOT40 +
            cloud + cloud ^ 2 +
            ", add_terms(c(
  str_c("bin", 0:39),
  str_c("surface_", 1:9)
)), "|
            x_y") %>%
  as.formula()

fml_grid_year <- str_c("GOSIF_sum ~
            AOD + AOD ^ 2 +
            AOT40 +
            cloud + cloud ^ 2 +
            ", add_terms(c(
  str_c("bin", 0:39),
  str_c("surface_", 1:9)
)), "|
            x_y + year") %>%
  as.formula()

fml_county_year <- str_c("GOSIF_sum ~
            AOD + AOD ^ 2 +
            AOT40 +
            cloud + cloud ^ 2 +
            ", add_terms(c(
  str_c("bin", 0:39),
  str_c("surface_", 1:9)
)), "|
            county + year") %>%
  as.formula()

fml_county_year_trend <- str_c("GOSIF_sum ~
            AOD + AOD ^ 2 +
            AOT40 +
            cloud + cloud ^ 2 +
            ", add_terms(c(
  str_c("bin", 0:39),
  str_c("surface_", 1:9)
)), "|
            x_y + county[year]") %>%
  as.formula()

# climate -----------------------------------------------------------------

fml_tmax_step <- str_c("GOSIF_sum ~
            AOD + AOD ^ 2 +
            AOT40 +
            cloud + cloud ^ 2 +
            ", add_terms(c(
  str_c("step", 1:10),
  str_c("surface_", 1:9)
)), "|
            x_y[year]") %>%
  as.formula()

fml_tmax_poly3 <- str_c("GOSIF_sum ~
            AOD + AOD ^ 2 +
            AOT40 +
            cloud + cloud ^ 2 +
            ", add_terms(str_c("surface_", 1:9)), " +
            maxtmp + maxtmp ^ 2 + maxtmp ^ 3 |
            x_y[year]") %>%
  as.formula()

fml_tmax_poly5 <- str_c("GOSIF_sum ~
            AOD + AOD ^ 2 +
            AOT40 +
            cloud + cloud ^ 2 +
            ", add_terms(str_c("surface_", 1:9)), " +
            maxtmp + maxtmp ^ 2 + maxtmp ^ 3 + maxtmp ^ 4 + maxtmp ^ 5 |
            x_y[year]") %>%
  as.formula()

fml_tmax_spline3 <- str_c("GOSIF_sum ~
            AOD + AOD ^ 2 +
            AOT40 +
            cloud + cloud ^ 2 +
            ", add_terms(str_c("surface_", 1:9)), " +
            maxtmp_ns3_1 + maxtmp_ns3_2 + maxtmp_ns3_3 |
            x_y[year]") %>%
  as.formula()

fml_tmax_spline7 <- str_c("GOSIF_sum ~
            AOD + AOD ^ 2 +
            AOT40 +
            cloud + cloud ^ 2 +
            ", add_terms(str_c("surface_", 1:9)), " +
            ", add_terms(str_c("maxtmp_ns7_", 1:7)), " |
            x_y[year]") %>%
  as.formula()

fml_prep_poly3 <- str_c("GOSIF_sum ~
            AOD + AOD ^ 2 +
            AOT40 +
            cloud + cloud ^ 2 +
            prep + prep ^2 + prep ^ 3 +
            ", add_terms(str_c("bin", 0:39)), "|
            x_y[year]") %>%
  as.formula()

fml_prep_spline3 <- str_c("GOSIF_sum ~
            AOD + AOD ^ 2 +
            AOT40 +
            cloud + cloud ^ 2 +
            prep_ns3_1 + prep_ns3_2 + prep_ns3_3 +
            ", add_terms(str_c("bin", 0:39)), "|
            x_y[year]") %>%
  as.formula()

fml_surface_poly3 <- str_c("GOSIF_sum ~
            AOD + AOD ^ 2 +
            AOT40 +
            cloud + cloud ^ 2 +
            surface + surface ^2 + surface ^ 3 +
            ", add_terms(str_c("bin", 0:39)), "|
            x_y[year]") %>%
  as.formula()

fml_surface_spline3 <- str_c("GOSIF_sum ~
            AOD + AOD ^ 2 +
            AOT40 +
            cloud + cloud ^ 2 +
            surface_ns3_1 + surface_ns3_2 + surface_ns3_3 +
            ", add_terms(str_c("bin", 0:39)), "|
            x_y[year]") %>%
  as.formula()


fml_moisture_root <- str_c("GOSIF_sum ~
            AOD + AOD ^ 2 +
            AOT40 +
            cloud + cloud ^ 2 +
            ", add_terms(c(
  str_c("bin", 0:39),
  str_c("root_", 1:9)
)), "|
            x_y[year]") %>%
  as.formula()


fml_root_poly3 <- str_c("GOSIF_sum ~
            AOD + AOD ^ 2 +
            AOT40 +
            cloud + cloud ^ 2 +
            root + root ^2 + root ^ 3 +
            ", add_terms(str_c("bin", 0:39)), "|
            x_y[year]") %>%
  as.formula()

fml_root_spline3 <- str_c("GOSIF_sum ~
            AOD + AOD ^ 2 +
            AOT40 +
            cloud + cloud ^ 2 +
            root_ns3_1 + root_ns3_2 + root_ns3_3 +
            ", add_terms(str_c("bin", 0:39)), "|
            x_y[year]") %>%
  as.formula()


fml_cloud_poly3 <- str_c("GOSIF_sum ~
            AOD + AOD ^ 2 +
            AOT40 +
            cloud + cloud ^ 2 + cloud ^ 3 +
            ", add_terms(c(
  str_c("bin", 0:39),
  str_c("surface_", 1:9)
)), "|
            x_y[year]") %>%
  as.formula()

fml_cloud_poly5 <- str_c("GOSIF_sum ~
            AOD + AOD ^ 2 +
            AOT40 +
            cloud + cloud ^ 2 + cloud ^ 3 + cloud ^ 4 + cloud ^ 5 +
            ", add_terms(c(
  str_c("bin", 0:39),
  str_c("surface_", 1:9)
)), "|
            x_y[year]") %>%
  as.formula()

fml_cloud_spline3 <- str_c("GOSIF_sum ~
            AOD + AOD ^ 2 +
            AOT40 +
            cloud_ns3_1 + cloud_ns3_2 + cloud_ns3_3 +
            ", add_terms(c(
  str_c("bin", 0:39),
  str_c("surface_", 1:9)
)), "|
            x_y[year]") %>%
  as.formula()

fml_cloud_spline7 <- str_c("GOSIF_sum ~
            AOD + AOD ^ 2 +
            AOT40 +
            ", add_terms(str_c("cloud_ns7_", 1:7)), " +
            ", add_terms(c(
  str_c("bin", 0:39),
  str_c("surface_", 1:9)
)), "|
            x_y[year]") %>%
  as.formula()

fml_no_temp <- str_c("GOSIF_sum ~
            AOD + AOD ^ 2 +
            AOT40 +
            cloud + cloud ^ 2 +
            ", add_terms(str_c("surface_", 1:9)), "|
            x_y[year]") %>%
  as.formula()

fml_no_moisture <- str_c("GOSIF_sum ~
            AOD + AOD ^ 2 +
            AOT40 +
            cloud + cloud ^ 2 +
            ", add_terms(str_c("bin", 0:39)), "|
            x_y[year]") %>%
  as.formula()

fml_no_cloud <- str_c("GOSIF_sum ~
            AOD + AOD ^ 2 +
            AOT40 +
            ", add_terms(c(
  str_c("bin", 0:39),
  str_c("surface_", 1:9)
)), "|
            x_y[year]") %>%
  as.formula()

fml_no_climate <- str_c("GOSIF_sum ~
            AOD + AOD ^ 2 +
            AOT40 |
            x_y[year]") %>%
  as.formula()
