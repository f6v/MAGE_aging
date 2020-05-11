library(tidyverse)
library(futile.logger)
library(yaml)
source('./utils.R')

config_file_name <- commandArgs(trailingOnly = TRUE)[1]
config <- yaml.load_file(config_file_name)

if(!dir.exists(config$results_path)) dir.create(config$results_path)

duplicated_subjects <- load_duplicated_subjects(config$duplicated_subjects_path)

flog.logger("prepare_counts", INFO, appender = appender.file("prepare_counts.log"))

counts_for_chromosome <- function(chr) {
  mage_results <- load_mage_results(config$mage_results_path, chr)
  
  load_counts_with_phenotypes(config$counts_path, duplicated_subjects, chr) %>%
    filter(position %in% mage_results$locus)
}

phenotypes <- load_phenotypes(config$phenotypes_path)
all_chromosomes <- 1:22

flog.info('Started', name="prepare_counts")

all_chromosomes %>%
  map(counts_for_chromosome) %>%
  bind_rows() %>%
  write_csv(paste(config$results_path, "counts_with_mage_loci.csv", sep = ""))

flog.info('Stopped', name="prepare_counts")
