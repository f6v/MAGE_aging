library(data.table)

pA_filt <- 0.15
SE_filt <- 0.0075

load_counts <- function(counts_paths, duplicated_subjects, chr) {
  counts_for_chr_path <- paste(counts_paths, "counts_SNP_chr", chr, ".txt", sep = "")
  counts_col_names <- c("chromosome", "position", "ref_alleles", "variationtype", "dbSNP_ref", "gene", "A", "T", "C", "G", "sample")
  counts <- fread(counts_for_chr_path, col.names = counts_col_names)
  counts <- counts[variationtype == "SNV"]
  counts <- counts[!grepl("\\.", ref_alleles)]
  counts <- counts[(A + T + C + G) != 0]
  counts <- counts[!(sample %in% duplicated_subjects)]
  setkey(counts, position)
  
  return(counts)
}

load_counts_with_phenotypes <- function(counts_paths, duplicated_subjects, chr) {
  return(
    load_counts(counts_paths, duplicated_subjects, chr) %>%
      mutate(subject_id = stringr::str_match(sample, "[^-]*-[^-]*")) %>%
      inner_join(phenotypes, by = "subject_id") %>%
      filter(variationtype == "SNV")
  )
}

load_phenotypes <- function(phenotypes_path) {
  return(
    read_tsv(file = phenotypes_path, skip = 10) %>%
      select(SUBJID, AGE, SEX) %>%
      rename(subject_id = SUBJID, age = AGE, sex = SEX)
  )
}

load_mage_results <- function(mage_results_path, chr) {
  results <- read_csv(paste(mage_results_path, "mage_processsed_", chr, ".csv", sep = "")) %>%
    filter(seq_error < SE_filt & ref_allele_freq > pA_filt & ref_allele_freq < (1 - pA_filt))

  inbr_coeff_chr <- median(results$inbr_coeff)

  return(
    results %>%
      mutate(
        frac_heteroz = 2 * ref_allele_freq * (1 - ref_allele_freq) * (1 - inbr_coeff_chr),
        chr = chr
      )
  )
}

seqem_path_for_process <- function(base_path) {
  return(paste(base_path, "pid_", Sys.getpid(), "/", sep = ""))
}

make_plot <- function(counts, mage_results, current_locus) {
  samples_for_analysis <- get_heteroz_samples(counts, mage_results, current_locus)

  return(
    samples_for_analysis %>%
      ggplot(aes(x = age, y = ratio)) +
      geom_point() +
      theme_minimal() +
      geom_smooth(method='lm', formula= y ~ x)
  )
}

get_heteroz_samples <- function(counts, mage_results, current_locus) {
  locus_data <- mage_results %>% filter(locus == current_locus)
  samples_for_locus <- counts %>% filter(position == current_locus)
  ref_counts <- samples_for_locus[, locus_data$ref_allele]
  var_counts <- samples_for_locus[, locus_data$var_allele]
  samples_for_locus$ratio <- ifelse(ref_counts > var_counts, var_counts / ref_counts, ref_counts / var_counts)
  num_heteroz_samples <- round(dim(samples_for_locus)[1] * locus_data$frac_heteroz)
  
  return(
    samples_for_locus %>%
      arrange(desc(ratio)) %>%
      head(num_heteroz_samples)
  )
}

load_duplicated_subjects <- function(subjects_path) {
  return(
    read_csv(subjects_path, col_names = c("subject_id"))$subject_id
  )
}



