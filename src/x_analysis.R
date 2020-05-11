library(tidyverse)
library(futile.logger)
library(doParallel)
library(yaml)
source('./utils.R')
source('./model_utils.R')

config_file_name <- commandArgs(trailingOnly = TRUE)[1]
config <- yaml.load_file(config_file_name)

if(!dir.exists(config$results_path)) dir.create(config$results_path)

flog.logger("x_analysis", INFO, appender = appender.file("x_analysis.log"))

inactive_genes <- read_delim(config$gene_status_path, delim = ";") %>%
  filter(xci_status == "inactive") %>%
  select(gene) %>%
  pull()

duplicated_loci = read_csv(config$duplicated_loci_path)$position
duplicated_subjects <- load_duplicated_subjects(config$duplicated_subjects_path)
phenotypes <- load_phenotypes(config$phenotypes_path) %>%
  filter(sex == 2) # female subjects
all_counts <- load_counts_with_phenotypes(config$counts_path, duplicated_subjects, "X") %>%
  filter(!(position %in% duplicated_loci), gene %in% inactive_genes)
all_mage_results <- load_mage_results(config$mage_results_path, "X")
all_subject_genotypes <- read_csv(paste(config$subject_genotypes_path, "genotypes_chrX.csv", sep = "")) %>%
  gather("subject_id", "posterior_prob", -locus) %>%
  filter(!is.na(posterior_prob), !(locus %in% duplicated_loci)) %>%
  mutate(short_subject_id = str_extract(subject_id, "^\\w+-\\w+"))
subject_analysis_results <- read_csv(config$subject_analysis_results)

subjects_with_genotypes <- unique(all_subject_genotypes$subject_id)

cl <- makeCluster(config$paralellism)
registerDoParallel(cl)

analysis_results <- foreach(current_subject = iter(subjects_with_genotypes), .packages = c("tidyverse", "data.table"), .verbose = F) %dopar% {
  #flog.info('Started processing subject: %s', current_subject, name="subject_analysis")
  
  subject_data <- all_subject_genotypes %>%
    filter(subject_id == current_subject) %>%
    inner_join(all_mage_results, by = "locus")
  
  expected_heteroz <- subject_data %>%
    summarize(mean(2 * ref_allele_freq * (1 - ref_allele_freq))) %>%
    pull()

  f_subject <- subject_analysis_results %>%
    filter(subject_id == current_subject) %>%
    select(f_subject) %>%
    pull()
  
  counts_for_subject <- all_counts %>%
    filter(sample == current_subject) %>%
    inner_join(subject_data, by = c("position" = "locus"))
  
  if(dim(counts_for_subject)[1] == 0) {
    return(NULL)
  }
  
  counts_with_ratios <- counts_for_subject %>%
    pivot_longer(cols = A:G) %>%
    group_by(position) %>%
    filter(name == ref_allele | name == var_allele) %>%
    mutate(new_count_cols = c('ref_count', 'var_count')) %>%
    ungroup %>%
    select(position, value, new_count_cols) %>%
    pivot_wider(names_from = new_count_cols, values_from = value) %>%
    left_join(counts_for_subject, .) %>%
    mutate(ratio = ifelse(ref_count > var_count, var_count / ref_count, ref_count / var_count))
  
  n_hetroz_loci <- round(dim(counts_with_ratios)[1] * expected_heteroz * (f_subject))
  
  subject_median_ratio <- counts_with_ratios %>%
    arrange(desc(ratio)) %>%
    head(n_hetroz_loci) %>%
    summarise(median(ratio)) %>%
    pull()
  
  #flog.info('Finished processing subject: %s', current_subject, name="subject_analysis")
  
  return(
    list(
      "subject_id" = current_subject,
      "expected_heteroz" = expected_heteroz,
      "f_subject" = f_subject,
      "n_hetroz_loci" = n_hetroz_loci,
      "subject_median_ratio" = subject_median_ratio
    )
  )
}

stopCluster(cl)

analysis_results %>%
  bind_rows() %>%
  write_csv(paste(config$results_path, "X_analysis.csv", sep = ""))



