library(yaml)
library(tidyverse)
source('./utils.R')
source('./model_utils.R')

config_file_name <- commandArgs(trailingOnly = TRUE)[1]
config <- yaml.load_file(config_file_name)

if(!dir.exists(config$results_path)) dir.create(config$results_path)

phenotypes <- load_phenotypes(config$phenotypes_path)
library_sizes = read_csv(config$library_sizes_path)
duplicated_subjects <- load_duplicated_subjects(config$duplicated_subjects_path)

run_expression_analysis<- function(chr, locus, gene) {
  
  counts <- load_counts_with_phenotypes(config$counts_path, duplicated_subjects, chr)
  mage_results <- load_mage_results(config$mage_results_path, chr)
  
  samples_for_gene <- get_heteroz_samples(counts, mage_results, locus) %>%
    inner_join(library_sizes, by = "sample") %>%
    mutate(total_count = A + T + C + G) %>%
    mutate(norm_count = log((total_count / num_reads) * mean(num_reads)))
  
  fit <- lm(norm_count ~ age, data = samples_for_gene)
  p_value <- lm_p_value(fit)
  r_sq <- summary(fit)$r.squared
  slope <- fit$coefficients[[2]]
  n_samples<-nrow(samples_for_gene)
  
  plot <- samples_for_gene %>%
    ggplot(aes(x = age, y = norm_count)) +
    geom_point() +
    theme_minimal() +
    geom_smooth(method='lm', formula = y ~ x)
  
  plot_path <- paste(config$results_path, "plot_", chr, "_", locus, ".png", sep = "")
  ggsave(plot_path, plot, "png")
  
  return(list("locus" = locus, "chr" = chr, "gene" = gene, "p_value" = p_value, "r_sq" = r_sq, "slope" = slope, "n_samples" = n_samples))
}

all_stat_analysis_results <- list.files(config$stat_analysis_results_path, pattern = "\\.csv$", full.names = T) %>%
  map(~read_csv(.x)) %>%
  bind_rows()

significant_genes <- all_stat_analysis_results %>%
  filter(p_adj < 0.1, !is.na(gene)) %>%
  select(gene) %>%
  pull()

all_stat_analysis_results %>%
  filter(gene %in% significant_genes) %>%
  select(chr, locus, gene) %>%
  pmap(run_expression_analysis) %>%
  bind_rows() %>%
  write_csv(paste(config$results_path, "genes_expression", ".csv", sep = ""))
