## code to prepare `KIDNEY` dataset goes here
library(tidyverse)
dt <- read_csv("data-raw/kidney_data.csv")
KIDNEY <- dt %>%
  mutate(
    time = rep(c("t1", "t2"), times = nrow(.)/2),
    delta.name = rep(c("delta1", "delta2"), times = nrow(.)/2)
  ) %>%
  select(-c(age, disease, frailty)) %>%
  pivot_wider(
    names_from = c(time, delta.name),
    values_from = c(t, delta)
  ) %>%
  magrittr::set_colnames(
    c("id", "sex", "t1", "t2", "delta1", "delta2")
  ) %>%
  mutate(
    sex = if_else(sex == 2, 0, sex)
  ) %>%
  identity()

usethis::use_data(KIDNEY, overwrite = TRUE)
