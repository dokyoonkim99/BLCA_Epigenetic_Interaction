library(data.table, lib = "/cbica/home/manu148/R/x86_64-redhat-linux-gnu-library/3.5")
library(zoo, lib = "/cbica/home/manu148/R/x86_64-redhat-linux-gnu-library/3.5")
library(lmtest, lib = "/cbica/home/manu148/R/x86_64-redhat-linux-gnu-library/3.5")

args = commandArgs(trailingOnly=TRUE)
ouput_dir = args[1]
outlier_type = args[2]
outlier_threshold = as.double(args[3])
sge_num = as.integer(args[4]) - 1
sge_total = as.double(args[5])

map = data.frame(fread("../generate_combinations/transcript_methylation_mirna_combinations_tcga_xena", header = F), stringsAsFactors = F)

isoform_data = data.frame(fread("../../data/qc_centered_data_xena/isoform_log_centered_data"), stringsAsFactors = F)
methylation_data = data.frame(fread("../../data/qc_centered_data_xena/methylation_centered_data"), stringsAsFactors = F)
mirnaseq_data = data.frame(fread("../../data/qc_centered_data_xena/mirna_log_centered_data"), stringsAsFactors = F)
clinical = data.frame(fread("../../data/qc_centered_data_xena/clinical_data"), stringsAsFactors = F)

#m = map[1,]
#isoform.name = "BHMT2_ENST00000255192"
#methylation.name = "cg01902605"
#mirna.name = "hsa-miR-106b-5p"

is_outlier = function(x, type, threshold) {
  if (type == "mean") {
    m = mean(x)
    sd = sd(x)
    return(x < (m - threshold * sd) | x > (m + threshold * sd))
  } else if (type == "IQR") {
    max = quantile(x, 0.75, na.rm=TRUE) + (IQR(x, na.rm=TRUE) * threshold )
    min = quantile(x, 0.25, na.rm=TRUE) - (IQR(x, na.rm=TRUE) * threshold )
    return(x < min | x > max)
  } else if (type == "None") {
    return(rep(F, length(x)))
  }
  return(NA)
}

calculate.pval = function(m) {
  m = as.vector(unlist(m))
  isoform.name = m[1]
  methylation.name = m[2]
  mirna.name = m[3]
  ret = paste(isoform.name, methylation.name, mirna.name, sep="\t")
  
  isoform_row = which(isoform_data$V1 == isoform.name)
  methylation_row = which(methylation_data$V1 == methylation.name)
  mirna_row = which(mirnaseq_data$V1 == mirna.name)
  if (length(isoform_row) == 0 || length(methylation_row) == 0 || length(mirna_row) == 0) {
    p = paste(ret, paste("Error Code 1", length(isoform_row) == 0, length(methylation_row) == 0, 
                         length(mirna_row) == 0, sep=" "), sep = "\t")
    cat(paste(p, '\n', sep=""))
    return(p)
  }
  
  isoform.expression = unlist(isoform_data[isoform_row,-1])
  methylation = unlist(methylation_data[methylation_row,-1])
  mirna = unlist(mirnaseq_data[mirna_row,-1])
  
  d = data.frame(isoform.expression, methylation, mirna, clinical$gender, clinical$age_at_initial_pathologic_diagnosis ,clinical$ajcc_pathologic_tumor_stage, clinical$histological_grade)
  d = d[!apply(d, 1, anyNA),]
  d = d[!is_outlier(d$isoform.expression, outlier_type, outlier_threshold) & !is_outlier(d$methylation, outlier_type, outlier_threshold) & !is_outlier(d$mirna, outlier_type, outlier_threshold), ]
  formula = ""
  if (length(unique(d$clinical.ajcc_pathologic_tumor_stage)) > 1) {
    formula = paste(formula, " + as.factor(clinical.ajcc_pathologic_tumor_stage)", sep = '')
  }
  if (length(unique(d$clinical.gender)) > 1) {
    formula = paste(formula, " + as.factor(clinical.gender)", sep = '')
  }
  if (length(unique(d$clinical.histological_grade)) > 1) {
    formula = paste(formula, " + as.factor(clinical.histological_grade)", sep = '')
  }
  
  full = lm(formula = as.formula(paste("isoform.expression ~ methylation + mirna + methylation*mirna + clinical.age_at_initial_pathologic_diagnosis", formula, sep = '')), data = d)
  reduced = lm(formula = as.formula(paste("isoform.expression ~ methylation + mirna + clinical.age_at_initial_pathologic_diagnosis", formula, sep = '')), data = d)
  output = lrtest(full, reduced)
  val = output$`Pr(>Chisq)`[2]
  coef_names = names(full$coefficients)
  beta_methylation = full$coefficients[[which(coef_names == "methylation")]]
  beta_mirna = full$coefficients[[which(coef_names == "mirna")]]
  beta_methy_mirna = full$coefficients[[which(coef_names == "methylation:mirna")]]
  
  mirna_lm = lm(formula = as.formula(paste("isoform.expression ~ mirna + clinical.age_at_initial_pathologic_diagnosis", formula, sep = '')), data = d)
  mirna_lm_su = summary(mirna_lm)
  mirna_pval = mirna_lm_su$coefficients[2,4]
  
  me_lm = lm(formula = as.formula(paste("isoform.expression ~ methylation + clinical.age_at_initial_pathologic_diagnosis", formula, sep = '')), data = d)
  me_lm_su = summary(me_lm)
  me_pval = me_lm_su$coefficients[2,4]
  
  corr_isoform_methyl = cor(d$isoform.expression, d$methylation, method = "spearman")
  corr_isoform_mirna = cor(d$isoform.expression, d$mirna, method = "spearman")
  corr_methyl_mirna = cor(d$methylation, d$mirna, method = "spearman")
  
  ret = paste(ret, nrow(d), val, me_pval, mirna_pval, beta_methylation, beta_mirna, beta_methy_mirna, 
              corr_isoform_methyl, corr_isoform_mirna, corr_methyl_mirna, sep = "\t")
  #p = paste("Output_tcga:", ret, '\n', sep="")
  #cat(p)
  return(ret)
}

p_start = sge_num / sge_total
p_end = (sge_num + 1) / sge_total
start = floor(nrow(map) * p_start) + 1
end = floor(nrow(map) * p_end)

pvals = apply(map[start:end,], 1, calculate.pval)
write.table(data.frame(pvals), file =paste(ouput_dir, sge_num, sep = .Platform$file.sep), quote = F, row.names = F, col.names = F)

