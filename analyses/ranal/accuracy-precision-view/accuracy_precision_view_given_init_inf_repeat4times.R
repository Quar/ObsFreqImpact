require(tidyverse)
require(scales)

plot.accuracy.precision = function(sampling.method=c("snapshot", "upperbound"), measure=c("median", "mode", "mean")) {

  source("util/configurations.R")
  source("util/preproc_accuracy_precision_given_init_inf.R")
  
  fpath = fpaths.incidence.summary.newpara.repeat4times[sampling.method]  

  ( preproc.accuracy.precision(fpath)
    %>% calc.measure.deviation.ratio
    %>% mutate(dataset = dataset %>% toupper %>% as_factor)
  ) -> df.obs.var.ratio
  
  measure_incidence_ratio = sym(paste0(measure, "_incidence_ratio"))
  
  (ggplot(df.obs.var.ratio, 
          aes(
              x=!!measure_incidence_ratio,
              y=iqr_incidence_ratio,
              ))
    + geom_point(aes(
      shape=sample_interval
      , color=sample_interval
      ),
      size = 1)
    + scale_color_manual(values = hue_pal()(6))
    + guides(
      color=guide_legend(title="Duty Cycle Interval")
      , shape=guide_legend(title="Duty Cycle Interval")
    )    
    + facet_grid(disease ~ dataset) 
    + geom_vline(xintercept=0, color="grey") 
    + geom_hline(yintercept=0, color="grey") 
    + theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            aspect.ratio=1,
            axis.text=element_text(size=9),
            strip.text=element_text(size=12)
            ) 
    + scale_x_continuous(limits=c(-1,1))
    + scale_y_continuous(limits=c(-1,1))
    + xlab(paste(stringr::str_to_title(measure) ,"Attack Rate Deviation (Accuracy)"))
    + ylab("IQR Attack Rate Deviation (Precision)")
    + scale_color_discrete(name="Sampling Interval")
  ) -> g
    
  return(g)

}



if (sys.nframe() <= 1) {
  
  cat("==> draw plots ...\n")
  
  scenarios = expand.grid(
    sampling.method = c("snapshot", "upperbound"),
    measure =  c("median")
  )
  
  for (scn.idx in 1:nrow(scenarios)) {
    
    scn = scenarios[scn.idx,]
    
    g = plot.accuracy.precision(sampling.method=scn$sampling.method, measure=scn$measure)

    fname = paste0("accuracy-precision-view/accuracy_precision_view_repeat4times_given_init_inf_", scn$sampling.method, "_", scn$measure, ".pdf")
    
    cat("--> write to ", fname, "\n")
    
    ggsave(fname, g, width=12, height=16)
  }
  
}



