library(tidyverse)
library(futile.logger)
library(doParallel)
library(yaml)
source('./utils.R')
source('./model_utils.R')

config_file_name <- commandArgs(trailingOnly = TRUE)[1]
config <- yaml.load_file(config_file_name)

if(!dir.exists(config$results_path)) dir.create(config$results_path)

phenotypes <- load_phenotypes(config$phenotypes_path)
if (config$female_only) {
  phenotypes <- phenotypes %>%
    filter(sex == 2)
}

duplicated_subjects <- load_duplicated_subjects(config$duplicated_subjects_path)

flog.logger("stat_analysis", INFO, appender = appender.file('stat_analysis.log'))

for(chr in config$chromosomes) {
  flog.info('Started processing for chromosome %s', chr, name="stat_analysis")
  
  counts <- load_counts_with_phenotypes(config$counts_path, duplicated_subjects, chr)
  mage_results <- load_mage_results(config$mage_results_path, chr)
  
  cl <- makeCluster(config$paralellism)
  registerDoParallel(cl)
  
  analysis_results <- foreach(current_locus = iter(mage_results$locus), .packages = c("tidyverse", "data.table")) %dopar% {
    samples_for_analysis <- get_heteroz_samples(counts, mage_results, current_locus)
    
    return(
      fit_model_for_locus(current_locus, samples_for_analysis)
    )
  }
  
  stopCluster(cl)
  
  analysis_results <- bind_rows(lapply(analysis_results, as.data.frame.list))
  analysis_results$p_adj <- p.adjust(analysis_results$p_value, method = "BH")
  
  analysis_results_sig <- analysis_results %>%
    filter(p_adj < 0.1)
  
  for(current_locus in analysis_results_sig$locus) {
    plot <- make_plot(counts, mage_results, current_locus)
    plot_path <- paste(config$plots_path, "plot_", chr, "_", current_locus, ".png", sep = "")
    ggsave(plot_path, plot, "png")
  }

  write_csv(analysis_results, path = paste(config$results_path, "stat_analysis_", chr, ".csv", sep = ""))
  
  flog.info('Saved output for chromosome %s', chr, name = "stat_analysis")
}
