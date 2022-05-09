require(dplyr)
require(stringr)

INCIDENCE_SUMMARY_OUTPUT_ROOT = '/path/to/incidence/summary/output/root'
INFECTION_RECORDS_OUTPUT_ROOT = '/path/to/infection/records/output/root'


INFECTION_RECORDS_TEMPLATE = c(
  '{INFECTION_RECORDS_OUTPUT_ROOT}',
  '{root.name[sampling.method]}/',
  '/rssi_-80/infection_pairs/',
  'data_{data.collection}/',
  'pp_{dataset}_{data.collection}/',
  'result_{experiment.id}_samp{observation.interval}_only1.csv'
) %>% str_c(collapse="")


fpaths.incidence.summary.newpara.repeat4times = c(
  snapshot=str_glue('{INCIDENCE_SUMMARY_OUTPUT_ROOT}/incidence_summary_newpara_repeat4times_snapshot.csv'),
  upperbound=str_glue('{INCIDENCE_SUMMARY_OUTPUT_ROOT}/incidence_summary_newpara_repeat4times_upperbound.csv')
)


fpaths.infection.records.only1 = function(sampling.method, data.collection, dataset, observation.interval, experiment.id) {
  
  stopifnot(sampling.method %in% c('snapshot', 'upperbound'))
  
  stopifnot(data.collection %in% c('origin', 'repeat4times'))
  
  stopifnot(dataset %in% paste0('shed', c(1,2,7,8,9)))
  
  stopifnot(observation.interval %in% c(5, 10, 30, 60, 90, 180, 360))
  
  root.name=c(
    snapshot='pp_newpara_snapshot',
    upperbound='pp_newpara_upperbound'
  )
  
  fpath = str_glue(INFECTION_RECORDS_TEMPLATE)
  
  stopifnot(file.exists(fpath))
  
  fpath
}

tbl.diseases.experiment.id = tibble::tribble(
  ~experiment.id, ~disease, 
  "e1", "flu",
  "e2", "sars",
  "e3", "fifth",
  "e4", "pertussis",
  "e5", "chickenpox",
  "e6", "mers",
  "e7", "diphtheria",
  "e8", "measles",
  "e9", "covid19",
  "e10", "covid19alpha",
  "e11", "covid19beta",
  "e12", "covid19delta"
)

DISEASES = tbl.diseases.experiment.id$disease

get.experiment.id = function(disease.name) {
  experiment.id = tbl.diseases.experiment.id %>% filter(disease==disease.name) %>% pull("experiment.id")
  has.found.experiment.id = length(experiment.id) > 0
  stopifnot( has.found.experiment.id )
  experiment.id
}

get.disease.name = function(experiment.index) {
  disease.name = tbl.diseases.experiment.id %>% filter(experiment.id==experiment.index) %>% pull("disease")
  has.found.disease.name = length(disease.name) > 0
  stopifnot( has.found.disease.name )
  disease.name
}