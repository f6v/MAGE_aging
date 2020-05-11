library(tidyverse)
library(futile.logger)
library(doParallel)
library(yaml)
source('./utils.R')
source('./model_utils.R')

## SETUP ##
config_file_name <- commandArgs(trailingOnly = TRUE)[1]
config <- yaml.load_file(config_file_name)

if(!dir.exists(config$results_path)) dir.create(config$results_path)

flog.logger("subject_analysis", INFO, appender = appender.file("subject_analysis.log"))
## END SETUP ##

duplicated_loci = read_csv(config$duplicated_loci_path)$position

all_counts <- read_csv(config$counts_with_mage_loci_path) %>%
  filter(!(position %in% duplicated_loci))

all_chromosomes <- 1:22
all_mage_results <- all_chromosomes %>%
  map(~load_mage_results(config$mage_results_path, .x)) %>%
  bind_rows() %>%
  filter(!(locus %in% duplicated_loci))

all_subject_genotypes <- all_chromosomes %>%
  map(~read_csv(paste(config$subject_genotypes_path, "genotypes_chr", .x, ".csv", sep = ""))) %>%
  bind_rows() %>%
  gather("subject_id", "posterior_prob", -locus) %>%
  filter(!is.na(posterior_prob), !(locus %in% duplicated_loci)) %>%
  mutate(short_subject_id = str_extract(subject_id, "^\\w+-\\w+"))

subjects_with_genotypes <- unique(all_subject_genotypes$subject_id)
## END LOAD DATA ##

cl <- makeCluster(config$paralellism)
registerDoParallel(cl)

analysis_results <- foreach(current_subject = iter(subjects_with_genotypes),
                            .packages = c("tidyverse", "data.table"), .verbose = F) %dopar% {
  #flog.info('Started processing subject: %s', current_subject, name="subject_analysis")

  subject_data <- all_subject_genotypes %>%
    filter(subject_id == current_subject) %>%
    inner_join(all_mage_results, by = "locus")
  
  actual_heteroz <- subject_data %>%
    summarize(mean(posterior_prob)) %>%
    pull()
  expected_heteroz <- subject_data %>%
    summarize(mean(2 * ref_allele_freq * (1 - ref_allele_freq))) %>%
    pull()
  f_subject <- actual_heteroz / expected_heteroz
  
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
      "actual_heteroz" = actual_heteroz,
      "f_subject" = f_subject,
      "n_hetroz_loci" = n_hetroz_loci,
      "subject_median_ratio" = subject_median_ratio
    )
  )
}

stopCluster(cl)

analysis_results %>%
  bind_rows() %>%
  write_csv(paste(config$results_path, "subject_analysis.csv", sep = ""))