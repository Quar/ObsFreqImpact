require(tidyverse)
require(glue)

plot.kld.facet = function(ifile, ofile=NULL) {
  if (is.null(ofile)) {
    ofile = gsub(".csv", ".pdf", ifile)
  }
  
  cat(glue("read from {ifile}, output to {ofile} ..."), sep = "\n")
  df = read.csv(ifile) %>% 
    mutate(dataset = dataset %>% toupper %>% as_factor)
  cat(glue("maximum y-value is: {max(df$kld)}"))
  pp = (ggplot(df)
        + aes(x=obs_interval, y=kld, color=dataset)
        + geom_line()
        + facet_wrap('disease')
        + theme(
            axis.text.x = element_text(angle=45, size=12),
            axis.title = element_text(size=14),
            legend.title = element_text(size=10) #,
        ) 
        + ylim(0, 3)
        + labs(x='Duty Cycle Interval (unit: minute)'
             , y='KL-Divergence')
  )
  
  ggsave(ofile, pp, width=7, height=5)
}


downsample.method = c("snapshot", "upperbound")
data.collection = c("repeat4times")

ifiles = crossing(downsample.method, data.collection) %>%
  pmap_dfr(~setNames(glue("kld/data/kld_{..1}_{..2}.csv"), "filename"))


df.all = ifiles %>% pmap_dfr(~read.csv(.))

for (i in ifiles %>% as.matrix) {
  dir.name = dirname(i)
  file.name = basename(i)
  output.dirname = "./kld"
  output.filename = gsub(".csv", ".pdf", file.name)
  ofile = file.path(output.dirname, output.filename)
  plot.kld.facet(i, ofile)
}