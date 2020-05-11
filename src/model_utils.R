lm_p_value <- function(modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm'")
  
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1], f[2], f[3], lower.tail = F)
  attributes(p) <- NULL
  
  return(p)
}

fit_model_for_locus <- function(locus, data) {
  r_sq <- NULL
  r_sq_adj <- NULL
  f_stat <- NULL
  slope <- NULL
  p_value <- NULL
  
  n_samples <- dim(data)[1]
  
  if(n_samples >= 8) {
    fit <- lm(data$ratio ~ data$age)
    
    r_sq <- summary(fit)$r.squared
    r_sq_adj <- summary(fit)$adj.r.squared
    f_stat <- summary(fit)$fstatistic[[1]]
    slope <- fit$coefficients[[2]]
    p_value <- lm_p_value(fit)
  }
  
  chr = data[1, "chromosome"]
  dbsnp_ref = data[1, "dbSNP_ref"]
  gene = data[1, "gene"]
  
  return(c("chr" = chr, "locus" = locus, "r_sq" = r_sq, "r_sq_adj" = r_sq_adj, "f_stat" = f_stat, "slope" = slope, "p_value" = p_value,
           "n_heteroz" = n_samples, "dbSNP_ref" =  dbsnp_ref, "gene" = gene))
}