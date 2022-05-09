library(ggplot2)
library(tidyverse)
source("util/format_helper.R")

library(showtext)
# font_add_google("Arial", "sans")
font_add("Arial", 
         regular="Arial.ttf",
         bold="Arial Bold.ttf",
         italic="Arial Italic.ttf",
         bolditalic = "Arial Bold Italic.ttf"
         )

showtext_auto()


plot.tot_incidence_num = function(fpath, gpath) {
  df = read.csv(fpath, stringsAsFactors = T)
  
  pop.size = c(shed1=39, shed2=32, shed7=61, shed8=74, shed9=78)
  
  df.pop.size = data.frame(
    dataset = df$dataset %>% levels %>% factor,
    pop_size = pop.size %>% as.numeric
  )
  
  ##### plot #####
  
  df$sample_interval = factor(df$sample_interval, c(5,10,30,60,90,180,360), ordered = T)
  df = df %>% format_names_dataset #%>% format_names_disease
  df.pop.size = df.pop.size %>% format_names_dataset
  
  ggplot(df) +
    geom_violin(aes(x=sample_interval, y=num_incidence), scale="width") +
    stat_summary(aes(x=sample_interval, y=num_incidence), fun=mean, geom="point", shape=20, size=2, color="red", fill="red") +    
    geom_hline(aes(yintercept=pop_size),  color="blue", df.pop.size) +
    facet_grid(disease ~ dataset, scales = "free") +
    xlab("Duty Cycle Interval (unit: minute)") +
    ylab("Cumulative Cases") +
    theme_grey(base_size = 8) +
    theme(
      # strip.text=element_text(size=8),
      text = element_text(family = "Arial", size=8)
    )
  
  ggsave(gpath, width=7.25, height=8.7)
}


scenarios = list(
  repeat4times = list(
    fpath = fpaths.incidence.summary.newpara.repeat4times['upperbound'],
    # gpath = 'violin-plot-incidence-num/violin_incidence_num_newpara_upperbound_repeat4times.pdf'
    gpath = 'violin-plot-incidence-num/violin_incidence_num_newpara_upperbound_repeat4times-portrait.pdf'
  ),
  repeat4times.tiff = list(
    fpath = fpaths.incidence.summary.newpara.repeat4times['upperbound'],
    # gpath = 'violin-plot-incidence-num/violin_incidence_num_newpara_upperbound_repeat4times.png'
    gpath = 'violin-plot-incidence-num/violin_incidence_num_newpara_upperbound_repeat4times-portrait.tiff'
  ),
  repeat4times.eps = list(
    fpath = fpaths.incidence.summary.newpara.repeat4times['upperbound'],
    # gpath = 'violin-plot-incidence-num/violin_incidence_num_newpara_upperbound_repeat4times.png'
    gpath = 'violin-plot-incidence-num/violin_incidence_num_newpara_upperbound_repeat4times-portrait.eps'
  )
)


for (scn in scenarios) {
  cat(str(scn))
  plot.tot_incidence_num(scn$fpath, scn$gpath)  
}