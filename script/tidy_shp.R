source("script/loadPackages.R")
library(rvest)
library(fs)

# 爬取省代码列表
codedf <-
  read_html("http://www.mca.gov.cn/article/sj/xzqh/2020/20201201.html") %>%
  html_nodes(xpath = "/html/body/div/table") %>%
  html_table() %>%
  .[[1]] %>%
  filter(X1 == "" & X2 != "") %>%
  select(X2, X3) %>%
  set_names(c("code", "name")) %>%
  slice(-1)

API_pre <- "http://xzqh.mca.gov.cn/data/"
## province-level
province <- st_read(dsn = paste0(API_pre, "quanguo.json")) %>%
  `st_crs<-`(4326) %>%
  select(QUHUADAIMA) %>%
  set_names(c("province_code", "geometry")) %>%
  left_join(codedf, c("province_code" = "code")) %>%
  rename(province = name)

## city-level
city <- province %>%
  st_drop_geometry() %>%
  filter(province_code != "daodian") %>%
  distinct() %>%
  mutate(data = map(province_code, safely(function(acode) {
    X <- st_read(dsn = paste0(API_pre, acode, ".json"))
    st_crs(X) <- 4326

    X %>%
      filter(!FillColor == "") %>%
      filter(!NAME == "") %>%
      select(NAME, QUHUADAIMA) %>%
      set_names(c("city", "city_code", "geometry"))
  })))

city <- city %>%
  mutate(flag = map_lgl(data, ~ is.null(.x$error))) %>%
  filter(flag) %>%
  unnest_wider(data) %>%
  select(province_code, province, result) %>%
  unnest()

city <- city %>%
  filter(str_sub(province_code, end = 2) == str_sub(city_code, end = 2))

city <- st_as_sf(city, sf_column_name = "geometry")

## county-level
county <- st_read(dsn = paste0(API_pre, "xian_quanguo.json")) %>%
  `st_crs<-`(4326) %>%
  select(NAME, QUHUADAIMA) %>%
  set_names(c("county", "county_code", "geometry")) %>%
  mutate(code1 = str_sub(county_code, 1, 2)) %>%
  left_join(codedf %>%
    filter(str_detect(code, "0000")) %>%
    mutate(code1 = str_sub(code, 1, 2)) %>%
    rename(province = name, province_code = code)) %>%
  select(-code1)

## jiuduanxian
jiuduanxian <- st_read(dsn = paste0(API_pre, "quanguo_Line.geojson")) %>%
  `st_crs<-`(4326) %>%
  filter(QUHUADAIMA == "guojiexian") %>%
  filter(NAME != "jiantou") %>%
  head(14)

# ggplot() +
#   geom_sf(data = county) +
#   geom_sf(data = China_line)

# st_write(province, 'data/inputs/shp/province.shp',
#          layer_options = "ENCODING=UTF-8")
# st_write(county, 'data/inputs/shp/county.shp',
#          layer_options = "ENCODING=UTF-8")
# st_write(jiuduanxian, 'data/inputs/shp/jiuduanxian.shp',
#          layer_options = "ENCODING=UTF-8")

write_rds(province, "data/inputs/shp/province.rds")
write_rds(city, "data/inputs/shp/city.rds")
write_rds(county, "data/inputs/shp/county.rds")
write_rds(jiuduanxian, "data/inputs/shp/jiuduanxian.rds")
