library(tidyverse)
library(futile.logger)
library(data.table)
library(MAGE)
library(yaml)
library(doParallel)
source('./utils.R')

pA_estimate_filt <- 0.1
#cov_filt <- 4
nr_samples_filt <- 30

config_file_name <- commandArgs(trailingOnly = TRUE)[1]
config <- yaml.load_file(config_file_name)

cov_filt <- config$coverage_filter

duplicated_subjects <- load_duplicated_subjects(config$duplicated_subjects_path)

if(!dir.exists(config$results_path)) dir.create(config$results_path)

unlink(config$base_seqem_path, recursive = T)
dir.create(config$base_seqem_path)

flog.logger("mage_processing", INFO, appender = appender.file('mage_processing.log'))

for(chr in config$chromosomes) {
  flog.info('Started processing for chromosome %s', chr, name="mage_processing")
  
  counts <- load_counts(config$counts_paths, duplicated_subjects, chr)
  loci <- unique(counts$position)
  
  cl <- makeCluster(config$paralellism)
  registerDoParallel(cl)
  
  results_list <- foreach(locus = iter(loci), .packages = c("data.table", "MAGE")) %dopar% {
    #flog.info('Processing locus %s', locus, name="mage_processing")
    setDTthreads(2)
    samples_for_locus <- as.data.frame(counts[position == locus])
    # Add reference and variant alleles
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
    
    seqem_res <- MAGE::estimate_parameters(samples_filt$ref_count,
                                           samples_filt$var_count,
                                           seqem_path,
                                           ref_allele = samples_filt[1, "ref"],
                                           var_allele = samples_filt[1, "var"],
                                           position = locus,
                                           dbSNP = samples_filt[1, "dbSNP_ref"],
                                           chr = chr)

    return(
      list(
        "locus" = locus,
        "ref_allele_freq" = seqem_res$allelefreq,
        "seq_error" = seqem_res$SE,
        "inbr_coeff" = seqem_res$inbr,
        "ref_allele" = samples_filt[1, "ref"],
        "var_allele" = samples_filt[1, "var"]
      )
    )
  }
  stopCluster(cl)

  results_for_chr <- rbindlist(results_list)
  fwrite(results_for_chr, paste(config$results_path, "mage_processsed_", chr, ".csv", sep = ""))
  
  flog.info('Saved output for chromosome %s', chr, name="mage_processing")
}
