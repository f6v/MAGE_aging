library(tidyverse)
library(futile.logger)
library(data.table)
library(MAGE)
library(yaml)
library(doParallel)
source('./utils.R')
source('./seqem_utils.R')

pA_estimate_filt <- 0.1
cov_filt <- 4
nr_samples_filt <- 30

config_file_name <- commandArgs(trailingOnly = TRUE)[1]
config <- yaml.load_file(config_file_name)

if(!dir.exists(config$results_path)) dir.create(config$results_path)

unlink(config$base_seqem_path, recursive = T)
dir.create(config$base_seqem_path)

duplicated_subjects <- load_duplicated_subjects(config$duplicated_subjects_path)

flog.logger("export_subject_genotypes", INFO, appender = appender.file('export_subject_genotypes.log'))

for(chr in config$chromosomes) {
  flog.info('Started processing for chromosome %s', chr, name="export_subject_genotypes")
  
  mage_results <- load_mage_results(config$mage_results_path, chr)
  counts <- load_counts(config$counts_paths, duplicated_subjects, chr)[position %in% mage_results$locus]
  loci <- unique(counts$position)
  
  cl <- makeCluster(config$paralellism)
  registerDoParallel(cl)
  
  results_list <- foreach(locus = iter(loci), .packages = c("data.table", "MAGE")) %dopar% {
    # flog.info('Started processing for locus %s on chromosome %s', locus, chr, name="export_subject_genotypes")
    
    samples_for_locus <- as.data.frame(counts[position == locus])
    samples_std <- MAGE::standard_alleles(samples_for_locus)
  
    samples_filt <- MAGE::prior_filter(samples_std,
                                       prior_allelefreq = pA_estimate_filt,
                                       coverage_filter = cov_filt,
                                       samples_filter = nr_samples_filt)
    if(is.null(samples_filt)) {
      return(NULL)
    }
    
    seqem_path <- seqem_path_for_process(config$base_seqem_path)
    if(!dir.exists(seqem_path)) dir.create(seqem_path)
    
    probs_heteroz <- estimate_genotypes(samples_filt$ref_count,
                                           samples_filt$var_count,
                                           seqem_path,
                                           ref_allele = samples_filt[1, "ref"],
                                           var_allele = samples_filt[1, "var"],
                                           position = locus,
                                           dbSNP = samples_filt[1, "dbSNP_ref"],
                                           chr = chr)
    
    probs_heteroz <- setNames(probs_heteroz, samples_filt$sample)
    probs_heteroz$locus = locus

    return(probs_heteroz)
  }

  stopCluster(cl)

  results_list %>%
    bind_rows() %>%
    write_csv(paste(config$results_path, "genotypes_chr", chr, ".csv", sep = ""))
  
  flog.info('Finished processing for chromosome %s', chr, name="export_subject_genotypes")
}






