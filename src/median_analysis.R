library(tidyverse)
library(futile.logger)
library(doParallel)
library(yaml)
source('./utils.R')
source('./model_utils.R')

config_file_name <- commandArgs(trailingOnly = TRUE)[1]
config <- yaml.load_file(config_file_name)

if(!dir.exists(config$results_path)) dir.create(config$results_path)

duplicated_subjects <- load_duplicated_subjects(config$duplicated_subjects_path)

flog.logger("stat_analysis", INFO, appender = appender.file('stat_analysis.log'))

for(chr in config$chromosomes) {
  flog.info('Started processing for chromosome %s', chr, name="stat_analysis")
  
  phenotypes <- load_phenotypes(config$phenotypes_path)
  
  if (chr == "X") {
    counts_path <- config$x_counts_path
    mage_results_path <- config$x_mage_results_path
    phenotypes <- phenotypes %>%
      filter(sex == 2)
  } else {
    counts_path <- config$counts_path
    mage_results_path <- config$mage_results_path
  }
  
  counts <- load_counts_with_phenotypes(counts_path, duplicated_subjects, chr)
  mage_results <- load_mage_results(mage_results_path, chr)
  
  cl <- makeCluster(config$paralellism)
  registerDoParallel(cl)
  
  analysis_results <- foreach(current_locus = iter(mage_results$locus), .packages = c("tidyverse", "data.table")) %dopar% {
    samples_for_analysis <- get_heteroz_samples(counts, mage_results, current_locus)
    
    return(
      list("locus" = current_locus, "median" = median(samples_for_analysis$ratio))
    )
  }
  
  stopCluster(cl)
  
  analysis_results <- bind_rows(lapply(analysis_results, as.data.frame.list))
  write_csv(analysis_results, path = paste(config$results_path, "median_analysis_", chr, ".csv", sep = ""))
  
  flog.info('Saved output for chromosome %s', chr, name = "stat_analysis")
}
