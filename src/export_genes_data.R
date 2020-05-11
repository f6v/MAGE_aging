library(yaml)
library(tidyverse)
source('./utils.R')
source('./model_utils.R')

load_significant_results <- function(file_name) {
  read_csv(file_name) %>% filter(p_adj < 0.1)
}

config_file_name <- commandArgs(trailingOnly = TRUE)[1]
config <- yaml.load_file(config_file_name)

if(!dir.exists(config$results_path)) dir.create(config$results_path)

phenotypes <- load_phenotypes(config$phenotypes_path)
library_sizes = read_csv(config$library_sizes_path)
duplicated_subjects <- load_duplicated_subjects(config$duplicated_subjects_path)

data_for_gene <- function(chr, gene) {
  samples_for_gene <- load_counts_with_phenotypes(config$counts_path, duplicated_subjects, chr) %>%
    filter(gene == gene)

  return(samples_for_gene)
}

list.files(config$stat_analysis_results_path, pattern = "\\.csv$", full.names = T) %>%
  map(load_significant_results) %>%
  bind_rows() %>%
  select(chr, gene) %>%
  filter(!is.na(gene)) %>%
  pmap(data_for_gene) %>%
  bind_rows() %>%
  write_csv(paste(config$results_path, "genes_for_analysis", ".csv", sep = ""))
