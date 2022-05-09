library(tidyverse)
library(ggpubr)
library(rstatix)
library(muStat)

source("util/configurations.R")

prentice_modified_friedman_test = function(incidence_summary_file) {

  df = read_csv(incidence_summary_file) %>%
    convert_as_factor(dataset, disease, init_node_idx, seed_number, sample_interval)
  
  df %>% head
  
  res = with(df %>% unite("block", 1:3),
     mu.friedman.test(num_incidence, sample_interval, block)
  )
  
  return(res)
}


fpaths = fpaths.incidence.summary.newpara.repeat4times


for (f in fpaths %>% names) {
  cat(f, "\n------------")
  print(prentice_modified_friedman_test(fpaths[[f]]) %>% t)
}


### uber group test

# N_MAX=3
N_MAX=Inf
df.all = fpaths %>% names %>% 
  map(~read_csv(fpaths[[.x]], n_max = N_MAX) %>% 
        mutate(
          sampling_method=str_split(.x, "\\.", simplify=T)[1],
          data_colletion=str_split(.x, "\\.", simplify=T)[2]
        ) %>%
        convert_as_factor(
          data_colletion,
          sampling_method,
          dataset,
          disease,
          init_node_idx
        )
  ) %>%
  { do.call(rbind, .) } %>%
  arrange(
    data_colletion,
    sampling_method,
    dataset,
    disease,
    init_node_idx
  )  %>%
  relocate(
    data_colletion, 
    sampling_method, 
    dataset, 
    disease, 
    init_node_idx, 
    sample_interval, 
    seed_number, 
    num_incidence
  )  
df.all

res.all = with(df.all %>% unite("block", c(-sample_interval, -seed_number, -num_incidence)),
           mu.friedman.test(num_incidence, sample_interval, block, correct=T)
)

res.all


res.all = with(df.all %>% 
                 unite("block", c(dataset, disease, init_node_idx)) %>%
                 unite("group", c(data_colletion, sampling_method, sample_interval))
                 ,
               mu.friedman.test(num_incidence, group, block, correct=T)
)

res.all


res.all = with(df.all %>% 
                 unite("group", c(
                   data_colletion, 
                   sampling_method, 
                   dataset, 
                   disease, 
                   init_node_idx, 
                   sample_interval, 
                   seed_number
                ))
               ,
               mu.wilcox.test(num_incidence, group, correct=T
                              , alternative="greater"
                              )
)

res.all