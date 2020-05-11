estimate_genotypes <- function (ref_counts, var_counts, wd_seqem, ref_allele = "A",
                                 var_allele = "T", position = 100, dbSNP = "rs100", chr = 1,
                                 HWE = 0, remove = TRUE, Debug = 0, MAX_iterations = 100,
                                 EM_threshold = 1e-08, Initial_error = 0.01, Fix_error = 0,
                                 Min_read_depth = 1)
{
  old <- setwd(tempdir())
  nr_samples <- length(ref_counts)
  text <- paste(paste(">chr", chr, sep = ""), dbSNP, position,
                ref_allele, var_allele, ref_counts + var_counts, var_counts,
                sep = " ")
  name_file <- paste(wd_seqem, "chr", chr, "_s", 1:nr_samples,
                     ".out", sep = "")
  out <- lapply(1:nr_samples, function(y) write(text[y], file = name_file[y]))
  f <- file(paste(wd_seqem, "imprinting_chr", chr, ".ctrl",
                  sep = ""), "w")
  writeLines(paste("Header_files: ", nr_samples, sep = ""),
             f)
  writeLines(name_file, f)
  writeLines(c(paste("HWE: ", HWE, sep = ""), paste("Debug: ",
                                                    Debug, sep = ""), paste("MAX_iterations: ", MAX_iterations,
                                                                            sep = ""), paste("EM_threshold: ", EM_threshold, sep = ""),
               paste("Initial_error: ", Initial_error, sep = ""), paste("Fix_error: ",
                                                                        Fix_error, sep = ""), paste("Min_read_depth: ", Min_read_depth,
                                                                                                    sep = ""), paste("Outfile: result_chr", chr, ".txt",
                                                                                                                     sep = "")), f)
  close(f)
  setwd(wd_seqem)
  system2("seqem", args = paste("imprinting_chr", chr, ".ctrl",
                                sep = ""), stdout = paste("resultcmd_chr", chr, ".txt",
                                                          sep = ""), stderr = paste("errorcmd_chr", chr, ".txt",
                                                                                    sep = ""))
  outputSeqEM <- readLines(paste("result_chr", chr, ".txt",
                                 sep = ""), n = -1)

  genotype_strings <- grep("genotype\\[[[:digit:]]+\\]= ", outputSeqEM, value = T)
  probs_heteroz <- read.table(text = gsub("[^0-9. ]", "", genotype_strings), fill = T)[, 3]
  
  if (remove == T) {
    junk <- dir(path = wd_seqem, pattern = paste("chr", chr,
                                                 sep = ""))
    removed <- file.remove(junk)
  }
  on.exit(setwd(old), add = TRUE)
  
  return(probs_heteroz)
}