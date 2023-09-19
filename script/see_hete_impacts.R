source("script/loadPackages.R")
source("script/loadFunctions.R")
source("script/loadFormulas.R")

a <- qread("data/impacts_ozone_AOT40.qs")
b <- qread("data/impacts_ozone_AOT40_hete.qs") %>%
  select(-c(tdata, idata))

a %>%
  slice(30) %>%
  pull() %>%
  pluck(1) %>%
  filter(peak_level == 60) %>%
  select(x, y, `50%`) %>%
  rast() %>%
  plot()

b %>%
  slice(30) %>%
  pull() %>%
  pluck(1) %>%
  filter(peak_level == 60) %>%
  select(x, y, `50%`) %>%
  rast() %>%
  plot()
