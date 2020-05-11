library(tidyverse)
library(yaml)
source('./utils.R')

config_file_name <- commandArgs(trailingOnly = TRUE)[1]
config <- yaml.load_file(config_file_name)

phenotypes <- load_phenotypes(config$phenotypes_path)
duplicated_subjects <- load_duplicated_subjects(config$duplicated_subjects_path)

extract_loci_from_files <- function(x) {
  chr = str_match(x, "plot_([0-9]+)_([0-9]+)")[1, 2]
  locus = str_match(x, "plot_([0-9]+)_([0-9]+)")[1, 3]
  
  return(list("chr" = chr, "locus" = locus))
}

extract_results <- function(chr_and_locus) {
  counts <- load_counts_with_phenotypes(config$counts_path, duplicated_subjects, chr_and_locus$chr)
  mage_results <- load_mage_results(config$mage_results_path, chr_and_locus$chr)
  
  counts_with_mage_results = counts %>%
    filter(position == chr_and_locus$locus) %>%
    inner_join(mage_results, by = c("position" = "locus"))
  
  return(counts_with_mage_results)
}

chrs_and_loci <- list.files(config$stat_analysis_results_path, pattern = "\\.png$") %>%
  map(extract_loci_from_files)

bind_rows(chrs_and_loci %>% map(extract_results)) %>%
  write_csv("sanity.csv")


