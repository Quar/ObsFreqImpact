library(tidyverse)
library(ggpubr)
library(rstatix)
library(broom)

source("util/configurations.R")

`%!in%` <- Negate(`%in%`)

fpaths = fpaths.incidence.summary.newpara.repeat4times

df.snapshot = fpaths['snapshot'] %>% read_csv() %>%
  convert_as_factor(dataset, disease, init_node_idx, seed_number, sample_interval)
df.snapshot %>% head


df.upperbound = fpaths['upperbound'] %>% read_csv() %>%
  convert_as_factor(dataset, disease, init_node_idx, seed_number, sample_interval)
df.upperbound %>% head

df.snapshot.lowR0 = df.snapshot %>% filter( disease %!in% c('chickenpox', 'measles', 'pertussis') )
df.snapshot.highR0 = df.snapshot %>% filter( disease %in% c('chickenpox', 'measles', 'pertussis') )
df.upperbound.lowR0 = df.upperbound %>% filter( disease %!in% c('chickenpox', 'measles', 'pertussis') )
df.upperbound.highR0 = df.upperbound %>% filter( disease %in% c('chickenpox', 'measles', 'pertussis') )

t.tests.over.obs.intervals = function(df) {
  t.test.output = function(obs.interval.a, obs.interval.b) {
    t.test(
      x = df %>% filter(sample_interval==obs.interval.a) %>% pull(num_incidence),
      y = df %>% filter(sample_interval==obs.interval.b) %>% pull(num_incidence)
    ) %>%
      tidy %>%
      mutate(
        t.test.x=obs.interval.a,
        t.test.y=obs.interval.b
      ) %>%
      select(t.test.x, t.test.y, everything())
  }
  
  target.obs.intervals = c(5, 10, 30, 60, 90, 180, 360)
  
  target.obs.intervals %>%
    combn(2) %>% as_tibble() %>%
    purrr::map_dfr(~t.test.output(.[1], .[2])) %>%
    list
}


t.test.results = df.snapshot.lowR0 %>%
  unite("blocks", 1:2) %>%
  group_by(blocks) %>%
  summarize(t.test.result = t.tests.over.obs.intervals(cur_data())) %>%
  unnest(t.test.result)

t.test.results %>% write_csv("welch-t-test-snapshot-lowR0.csv")


t.test.results = df.upperbound.lowR0 %>%
  unite("blocks", 1:2) %>%
  group_by(blocks) %>%
  summarize(t.test.result = t.tests.over.obs.intervals(cur_data())) %>%
  unnest(t.test.result)

t.test.results %>% write_csv("welch-t-test-upperbound-lowR0.csv")