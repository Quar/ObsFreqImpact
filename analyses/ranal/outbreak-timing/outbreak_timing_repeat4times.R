require(tidyverse)
source("util/format_helper.R")

source("util/configurations.R")
# import get.experiment.id(...)
# import fpaths.infection.records.only1(...)

cum.incidence.over.time = function(obs.interval=10, 
                                   data.collection=c('repeat4times', 'origin'), 
                                   dataset="shed1", 
                                   plot=F, 
                                   disease.name) {

  str_glue("--> Reading outbreak timing: disease={disease.name}, dataset={dataset}, obs.interval={obs.interval}") %>% print
  
  get.fpath = function(dataset.name, observation.interval, sampling.method=c("snapshot", "upperbound")) {
    
    experiment.id = get.experiment.id(disease.name)
    
    fpaths.infection.records.only1(
      sampling.method, 
      data.collection="repeat4times", 
      dataset.name, 
      observation.interval, 
      experiment.id
      )
  }

  # str_glue('\n--> reading {get.fpath(dataset, 5, "upperbound")} ...\n') %>% print
  df.ref = read.csv(get.fpath(dataset, 5, "upperbound"))
  
  # str_glue('\n--> reading {get.fpath(dataset, obs.interval, "upperbound")} ...\n') %>% print
  df.ub = read.csv(get.fpath(dataset, obs.interval, "upperbound"))
  
  # str_glue('\n--> reading {get.fpath(dataset, obs.interval, "snapshot")} ...\n') %>% print
  df.ss = read.csv(get.fpath(dataset, obs.interval, "snapshot"))

  df.ref.ecdf = ecdf(df.ref$timestampInMinutes)
  df.ub.ecdf = ecdf(df.ub$timestampInMinutes)
  df.ss.ecdf = ecdf(df.ss$timestampInMinutes)

  xs.ref = df.ref$timestampInMinutes %>% unique %>% sort
  xs.ub = df.ub$timestampInMinutes %>% unique %>% sort
  xs.ss = df.ss$timestampInMinutes %>% unique %>% sort

  ys.ref = df.ref.ecdf(xs.ref)
  ys.ub = df.ub.ecdf(xs.ub)
  ys.ss = df.ss.ecdf(xs.ss)
  
  ys.counts.ref = floor(ys.ref * nrow(df.ref) / nrow(df.ref)) 
  ys.counts.ub = floor(ys.ub * nrow(df.ub) / nrow(df.ref))
  ys.counts.ss = floor(ys.ss * nrow(df.ss) / nrow(df.ref))

  if(plot) {
    plot(xs.ub, ys.ub, col='blue', type='l' )
    lines(xs.ss, ys.ss, col='red')
    lines(xs.ref, ys.ref, col='black')
  }
  
  rbind(
    data.frame(
      time.in.minutes = xs.ref,
      cum.fraction.incidence = ys.ref,
      cum.counts.incidence = ys.counts.ref,
      label = "baseline",
      dataset = dataset,
      obs.interval = obs.interval
    ),
    data.frame(
      time.in.minutes = xs.ub,
      cum.fraction.incidence = ys.ub,
      cum.counts.incidence = ys.counts.ub,
      label = "upperbound",
      dataset = dataset,
      obs.interval = obs.interval      
    ),
    data.frame(
      time.in.minutes = xs.ss,
      cum.fraction.incidence = ys.ss,
      cum.counts.incidence = ys.counts.ss,
      label = "snapshot",
      dataset = dataset,
      obs.interval = obs.interval      
    )
  ) %>% as_tibble()
}


plot.cum.incidence.over.time.for.disease = function(disease.name, y.values=c('fraction', 'counts'), save.to.file=F) {
  
  stopifnot(disease.name %in% DISEASES)
  
  experiment.id = get.experiment.id(disease.name)
  
  datasets = paste0("shed", c(1,2,7,8,9))
  obs.intervals = c(10, 30, 60, 90, 180, 360)
  
  (
    expand.grid(
      d = datasets,
      o = obs.intervals
    )
    %>% rowwise
    %>% do (
      cum.incidence.over.time(
        obs.interval = .$o,
        data.collection = "repeat4times",
        dataset = .$d,
        plot = F,
        disease.name = disease.name
      )
    )
  ) -> df.summary 
  
  df.summary$disease = disease.name
  df.summary = df.summary %>% format_names_dataset %>% format_names_disease
  
  #if (echo) df.summary %>% head
  
  col.y.names = c(
    fraction='cum.fraction.incidence',
    counts='cum.counts.incidence'
  )
  
  col.y.name = sym(col.y.names[y.values])
  
  str_glue("\n--> Use {col.y.name} | ") %>% cat
  
  (
    ggplot(df.summary, aes(x=time.in.minutes, y = !!col.y.name, color = label))
    + geom_line(position=position_jitter(w=0.02, h=0))
    + facet_grid(obs.interval ~ dataset)
    + xlab("Incidence Occurrence Time (unit: Minute)")
    + ylab("ecdf(Incidence Occurrence Time)")
    + theme(axis.text.x = element_text(angle = 90, hjust = 1))
    + scale_color_discrete(name="Sampling\nMethod")
  ) -> g  
  
  if (save.to.file) {
    output.file.path = str_glue(
      'outbreak-timing/overall_outbreak_timing_repeat4times_',
      '{y.values}_{experiment.id}_{disease.name}',
      '.pdf'
    )
    
    str_glue("--> Write to {output.file.path}\n") %>% cat
    
    ggsave(
      output.file.path,
      g, 
      width=12, 
      height=10.8
    )
  }
  
  return(list(g, df.summary))
}


############# main ######

for ( disease.name in DISEASES ) {
  # plot.cum.incidence.over.time.for.disease(disease.name, y.values = "counts", save.to.file=T)
  plot.cum.incidence.over.time.for.disease(disease.name, y.values = "fraction", save.to.file=T)
}
