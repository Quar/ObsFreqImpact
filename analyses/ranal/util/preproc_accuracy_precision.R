require(tidyverse)

preproc.accuracy.precision = function(fpath) {
  
  cat("--> Reading from: ", fpath, " ...")
  
  df = read.csv(fpath)

  pop.size = c(shed1=39, shed2=32, shed7=61, shed8=74, shed9=78)

  df$sample_interval = factor(df$sample_interval, c(5,10,30,60,90,180,360), ordered = T)

  aggr.mode <- function(codes){
    which.max(tabulate(codes))
  }

  (df 
    %>% group_by(dataset, disease, sample_interval)
    %>% summarize(mean_incidence = mean(num_incidence)
                  , median_incidence = median(num_incidence)
                  , var_incidence = var(num_incidence)
                  , mode_incidence = aggr.mode(num_incidence)              
                  , iqr_incidence= diff(quantile(num_incidence, c(.25, .75))) )
  ) -> df.obs.var
  
  return(df.obs.var)
  
}


calc.measure.deviation.ratio = function(df.obs.var) {
  
  pop.size = c(shed1=39, shed2=32, shed7=61, shed8=74, shed9=78)
  
  (df.obs.var
   %>% group_by(dataset, disease)  
   %>% mutate(mean_incidence_diff = mean_incidence - mean_incidence[sample_interval==5],
              median_incidence_diff = median_incidence - median_incidence[sample_interval==5],
              mode_incidence_diff = mode_incidence - mode_incidence[sample_interval==5],                      
              var_incidence_diff = var_incidence - var_incidence[sample_interval==5],
              iqr_incidence_diff = iqr_incidence - iqr_incidence[sample_interval==5])
   %>% mutate(median_incidence_ratio = median_incidence_diff / pop.size[dataset],
              mode_incidence_ratio = mode_incidence_diff / pop.size[dataset],
              mean_incidence_ratio = mean_incidence_diff / pop.size[dataset],
              iqr_incidence_ratio = iqr_incidence_diff / pop.size[dataset])  
   %>% filter(sample_interval!=5)
  ) -> df.obs.var.ratio
  
  return(df.obs.var.ratio)
  
}