source("script/loadPackages.R")
source("script/loadFunctions.R")

order <- c(
  "North China", "Northeast China", "East China", "South China",
  "Central China", "Southwest China", "Northwest China"
)

province_plot <- st_read("https://geo.datav.aliyun.com/areas_v3/bound/100000_full.json") %>%
  mutate(region = fcase(
    name %in% c("辽宁省", "吉林省", "黑龙江省"), "Northeast China",
    name %in% c("上海市", "江苏省", "浙江省", "安徽省", "福建省", "江西省", "山东省", "台湾省"), "East China",
    name %in% c("北京市", "天津市", "河北省", "山西省", "内蒙古自治区"), "North China",
    name %in% c("河南省", "湖北省", "湖南省"), "Central China",
    name %in% c("广东省", "广西壮族自治区", "海南省", "香港特别行政区", "澳门特别行政区"), "South China",
    name %in% c("四川省", "贵州省", "云南省", "西藏自治区", "重庆市"), "Southwest China",
    name %in% c("陕西省", "甘肃省", "青海省", "宁夏回族自治区", "新疆维吾尔自治区"), "Northwest China"
  ))

region_plot <-
  province_plot %>%
  st_make_valid() %>%
  group_by(region) %>%
  summarise(geometry = st_union(geometry))

city <-
  province_plot %>%
  st_drop_geometry() %>%
  filter(adcode != "100000_JD") %>%
  transmute(
    adcode = fifelse(adcode %in% c("110000", "120000", "310000", "500000", "460000", "710000", "810000", "820000"), adcode, str_c(adcode, "_full")),
    data = map(adcode, ~ st_read(str_c("https://geo.datav.aliyun.com/areas_v3/bound/", .x, ".json")))
  ) %>%
  as_tibble()

city <-
  bind_rows(
    city %>% filter(str_detect(adcode, "full")) %>% unnest(),
    city %>% filter(!str_detect(adcode, "full")) %>% unnest()
  ) %>%
  st_as_sf(sf_column_name = "geometry")

data <- qread("data/tidied.qs", nthreads = qn) %>%
  mutate(
    crop = fifelse(str_detect(crop, "Rice"), "Rice", crop),
    AOT40 = fcase(
      crop == "Maize", AOT40_6,
      crop == "Rice", AOT40_7,
      crop == "Wheat", AOT40_7 - AOT40_1
    ),
    O3 = fcase(
      crop == "Maize",
      O3_6,
      crop == "Rice",
      O3_7,
      crop == "Wheat",
      (O3_7 * (MA - `GR&EM` + 1) - O3_1) / (MA - `GR&EM` + 0)
    )
  )

yearmean <- data %>%
  lazy_dt() %>%
  group_by(x, y, crop) %>%
  summarise(across(c(contains("GOSIF"), PM25, AOD, O3, AOT40, fraction), fmean), .groups = "drop") %>%
  pivot_wider(names_from = crop, values_from = c(contains("GOSIF"), PM25, AOD, O3, AOT40, fraction)) %>%
  as_tibble() %>%
  rast(type = "xyz", crs = "epsg:4326")

c("Maize", "Wheat", "Rice") %>%
  walk(
    function(acrop) {
      extracted <- exact_extract(yearmean[acrop][[1:7]], city, "weighted_mean", weights = yearmean[acrop][[8]]) %>%
        bind_cols(city, .) %>%
        set_names(~ str_remove(.x, "weighted_mean.")) %>%
        set_names(~ str_remove(.x, str_c("_", acrop))) %>%
        mutate(
          fraction = exact_extract(yearmean[acrop][[8]], city, "sum"),
          area = units::drop_units(st_area(.)) / 1e6,
          fraction = fraction / area * 100
        )

      patch <- pmap(
        list(
          c("GOSIF_sum", "fraction", "PM25", "AOD", "O3", "AOT40"),
          list(
            WrensBookshelf::WB_brewer("ThisMooseBelongsToMe", direction = -1),
            MetBrewer::met.brewer("VanGogh3"),
            NatParksPalettes::natparks.pals("Glacier", direction = -1),
            MetBrewer::met.brewer("OKeeffe2"),
            NatParksPalettes::natparks.pals("Arches2", direction = -1),
            WrensBookshelf::WB_brewer("WhatWellBuild", direction = -1)
          ),
          c("SIF (mW $m^{-2}$ $nm^{-1}$ $sr^{-1}$)", "Crop fraction (%)", "PM$_{2.5}$ ($\\mu$g m$^{-3}$)", "AOD", "MDA8 ($\\mu$g m$^{-3}$)", "AOT40 (ppb h)")
        ),
        function(avar, aple, aname) {
          extracted %>%
            ggplot(aes(fill = !!sym(avar))) +
            geom_sf(color = NA) +
            geom_sf(data = region_plot, fill = NA, size = 0.8) +
            scale_fill_gradientn(
              name = TeX(aname),
              colours = aple,
              limits = extracted %>% pull({{ avar }}) %>% quantile(c(0.1, 0.9), na.rm = T),
              oob = squish,
              na.value = "grey90",
              guide = guide_colorbar(
                title.position = "top", barwidth = 10,
                label.position = "bottom"
              )
            )
        }
      ) %>%
        wrap_plots(ncol = 2) &
        theme_ipsum_rc(
          panel_spacing = grid::unit(0, "lines"),
          plot_margin = margin(0, 0, 0, 0),
          grid = F, axis = F,
          base_size = 15,
          axis_title_size = 15,
          strip_text_size = 20,
        ) &
        theme(
          legend.position = c(0.2, 0.2),
          legend.direction = "horizontal",
          axis.text.x = element_blank(),
          axis.text.y = element_blank()
        )

      ggsave(str_c("figures/background_map_", acrop, ".pdf"), patch, width = 2, height = 3, scale = 7)
    }
  )
