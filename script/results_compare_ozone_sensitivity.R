source("script/loadPackages.R")
source("script/loadFunctions.R")
source('script/loadFormulas.R')

feng <-
  read_xlsx("data/inputs/Field_studies/FengNF2022.xlsx",
            sheet = "AOT40-RY", range = "A4:M1000"
  ) %>%
  janitor::remove_empty(c("rows", "cols")) %>%
  rename(AOT40 = `AOT40 (ppm h)`) %>% 
  select(AOT40, RY_inbred_rice, RY_wheat, RY_maize) %>% 
  `names<-`(c('AOT40', 'Rice', 'Wheat', 'Maize')) %>% 
  pivot_longer(-AOT40, names_to = 'crop', values_to = 'Feng et al.') %>% 
  mutate(`Feng et al.` = (`Feng et al.` - 1))

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

our = 
  map_dfr(model, tidy, .id = 'crop') %>% 
  filter(term == 'AOT40') %>% 
  select(crop, estimate) %>% 
  mutate(estimate = estimate * 1e3,
         estimate = map(estimate, ~ unique(feng$AOT40) * .x),
         AOT40 = map(estimate, ~ unique(feng$AOT40))) %>% 
  unnest() %>% 
  rename(Our = estimate)

p = 
  inner_join(our, feng) %>% 
  pivot_longer(-c(crop, AOT40)) %>% 
  ggplot(aes(x = AOT40, y = value, color = name)) +
  facet_wrap(~crop) +
  geom_line(size = 1.5) +
  scale_y_continuous(labels = label_percent()) +
  theme_half_open(15, font_family = "Roboto Condensed") +
  background_grid() +
  theme(legend.title = element_blank(),
        legend.position = c(0.2, 0.2)) +
  stat_poly_eq(use_label(c("eq")), formula = y ~ x - 1,
               label.x = 'right') +
  ylab('Relative yield') +
  xlab('AOT40 (ppm h)')

ggsave('figures/compare_sensitivity.pdf', p, width = 2, height = 1, scale = 4)


