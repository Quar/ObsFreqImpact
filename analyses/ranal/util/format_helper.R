require(stringr)

format_names_dataset = function(df) {
  dataset_names = df %>% pull(dataset) %>% unique %>% sort
  df %>% mutate(
    dataset = factor(
      dataset,
      levels = dataset_names,
      labels = dataset_names %>% str_to_upper
    )
  )
}

recode_disease_names = function(x) {
  x %>% 
    as.character %>%
    str_to_lower %>%
    recode(
      "flu" = "Flu",
      "norovius" = "Norovirus",
      "sars" = "SARS",
      'fifth' = "Fifth",
      "pertussis" = "Pertussis",
      "chickenpox" = "Chickenpox",
      "mers" = "MERS",
      "h7n9" = "H7N9",
      "diphtheria" = "Diphtheria",
      "meningococcal" = "Meningococcal",
      "measles" = "Measles",
      "covid19" = "COVID19",
      "covid19alpha" = "COVID19Alpha",
      "covid19beta" = "COVID19Beta",
      "covid19delta" = "COVID19Delta"
    )
}


format_names_disease = function(df) {
  disease_names = df %>% pull(disease) %>% unique %>% sort
  df %>% mutate(
    disease = factor(
      disease,
      levels = disease_names,
      labels = disease_names %>% recode_disease_names
    )
  )
}

