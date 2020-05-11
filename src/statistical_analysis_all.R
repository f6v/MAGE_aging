library(tidyverse)
library(futile.logger)
library(doParallel)
library(yaml)
source('./utils.R')
source('./model_utils.R')

config_file_name <- commandArgs(trailingOnly = TRUE)[1]
config <- yaml.load_file(config_file_name)

if(!dir.exists(config$results_path)) dir.create(config$results_path)
if(!dir.exists(config$plots_path)) dir.create(config$plots_path)

phenotypes <- load_phenotypes(config$phenotypes_path)
duplicated_subjects <- load_duplicated_subjects(config$duplicated_subjects_path)

flog.logger("stat_analysis", INFO, appender = appender.file('stat_analysis.log'))

significant_loci_results <- list.files(config$stat_analysis_results_path, pattern = "\\.csv$", full.names = T) %>%
  map(~read_csv(.x)) %>%
  bind_rows() %>%
  filter(p_adj < 0.1)

for(chr in significant_loci_results$chr) {
  flog.info('Started processing for chromosome %s', chr, name="stat_analysis")
  
  sig_genes_in_chr <- significant_loci_results %>%
    filter(chr == !!chr) %>%
    select(gene) %>%
    pull()

  counts <- load_counts_with_phenotypes(config$counts_path, duplicated_subjects, chr) %>%
    filter(gene %in% sig_genes_in_chr)
  mage_results <- load_mage_results(config$mage_results_path, chr)
  loci_in_gene = intersect(mage_results$locus, counts$position)
  
  cl <- makeCluster(config$paralellism)
  registerDoParallel(cl)
  
  analysis_results <- foreach(current_locus = iter(loci_in_gene), .packages = c("tidyverse", "data.table")) %dopar% {
    samples_for_analysis <- get_heteroz_samples(counts, mage_results, current_locus)
    
    return(
      fit_model_for_locus(current_locus, samples_for_analysis)
    )
  }
  
  stopCluster(cl)
  
  analysis_results <- bind_rows(lapply(analysis_results, as.data.frame.list))
  analysis_results$p_adj <- p.adjust(analysis_results$p_value, method = "BH")
  
  for(current_locus in analysis_results$locus) {
    plot <- make_plot(counts, mage_results, current_locus)
    plot_path <- paste(config$plots_path, "plot_", chr, "_", current_locus, ".png", sep = "")
    ggsave(plot_path, plot, "png")
  }

  write_csv(analysis_results, path = paste(config$results_path, "stat_analysis_", chr, ".csv", sep = ""))
  
  flog.info('Saved output for chromosome %s', chr, name = "stat_analysis")
}
