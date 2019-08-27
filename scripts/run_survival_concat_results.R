library(effects)
library(ggplot2)
library(survival)
library(stringr)
library(data.table)
library(GGally)

isoform_data = data.frame(fread("../data/qc_centered_data_xena/isoform_log_centered_data"), stringsAsFactors = F)
methylation_data = data.frame(fread("../data/qc_centered_data_xena/methylation_centered_data"), stringsAsFactors = F)
mirnaseq_data = data.frame(fread("../data/qc_centered_data_xena/mirna_log_centered_data"), stringsAsFactors = F)
clinical = data.frame(fread("../data/qc_centered_data_xena/clinical_data"), stringsAsFactors = F)
clinical$OS.time = clinical$OS.time /  30.417
clinical$surv = Surv(clinical$OS.time, clinical$OS)
outlier_type = "IQR" 
outlier_threshold = 3

# result file generated from run_lrt.R
map = read.csv(file = "xena_results_iqr3.tsv", sep = '\t', header = F, stringsAsFactors = F)
map$bonferroni = p.adjust(map$V5, method = "bonferroni", n = nrow(map))
map = map[map$bonferroni < 0.05,]
colnames(map) = c("Isoform", "methylation", "miRNA", "N", "pval_lrt", "pval_methylation", "pval_miRNA", "beta_methylation", 
                                  "beta_mirna", "beta_methy_mirna", "spearman_isoform_methyl", "spearman_isoform_mirna", "spearman_methyl_mirna", "bonferroni")
map = map[order(map$pval_lrt),]

#map = read.csv(file = 'results_outlier_iqr3_bonferroni_sig.csv', stringsAsFactors = F)
kp_pval = vector(length = nrow(map))
cox_pval = vector(length = nrow(map))
t_test = vector(length = nrow(map))
anova = vector(length = nrow(map))

isoform.name = "SGCD_ENST00000435422"
methylation.name = "cg19748027"
mirna.name = "hsa-miR-409-3p"

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

for (i in 1:nrow(map)) {
  m = as.vector(unlist(map[i,c(1:3)]))
  isoform.name = m[1]
  methylation.name = m[2]
  mirna.name = m[3]
  
  isoform_row = which(isoform_data$V1 == isoform.name)
  methylation_row = which(methylation_data$V1 == methylation.name)
  mirna_row = which(mirnaseq_data$V1 == mirna.name)
  
  isoform.expression = unlist(isoform_data[isoform_row,-1])
  methylation = unlist(methylation_data[methylation_row,-1])
  mirna = unlist(mirnaseq_data[mirna_row,-1])
  
  d = data.frame(isoform.expression, methylation, mirna, clinical$gender, clinical$age_at_initial_pathologic_diagnosis ,clinical$ajcc_pathologic_tumor_stage, clinical$histological_grade, clinical$surv)
  d = d[!apply(d, 1, anyNA),]
  d = d[!is_outlier(d$isoform.expression, outlier_type, outlier_threshold) & !is_outlier(d$methylation, outlier_type, outlier_threshold) & !is_outlier(d$mirna, outlier_type, outlier_threshold), ]
  
  methylation.breaks = c(quantile(d$methylation, probs = seq(0, 1, by = 0.333333333333333333333333333333333333333)))
  methylation.cut = cut(d$methylation, breaks=methylation.breaks, labels=c('1','2','3'), include.lowest=TRUE)
  
  mirna.breaks = c(quantile(d$mirna, probs = seq(0, 1, by = 0.333333333333333333333333333333333333333)))
  mirna.cut = cut(d$mirna, breaks=mirna.breaks, labels=c('1','2','3'), include.lowest=TRUE)
  
  quad = c(which(methylation.cut == 1 & mirna.cut == 1), which(methylation.cut == 3 & mirna.cut == 3))
  groups = c(rep(1, sum(methylation.cut == 1 & mirna.cut == 1)), rep(2, sum(methylation.cut == 3 & mirna.cut == 3)))
  
  surv_d = d[quad,]
  surv_d$groups = groups
  #surv_d = data.frame(surv[quad], groups)
  sdiff = survdiff(clinical.surv ~ groups, data = surv_d)
  fit  = survfit(clinical.surv ~ groups, data = surv_d)
  kp_pval[i] = 1 - pchisq(sdiff$chisq, length(sdiff$n) - 1)
  
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
  
  cox_reg = coxph(as.formula(paste("clinical.surv ~ groups + clinical.age_at_initial_pathologic_diagnosis + ", formula, sep = '')), data = surv_d)
  cox_summary = summary(cox_reg)
  cox_pval[i] = cox_summary$coefficients[1,5]
  
  t_test[i] = t.test(surv_d$isoform.expression ~ surv_d$groups)$p.value
  
  quad = c(which(methylation.cut == 1 & mirna.cut == 1), which(methylation.cut == 3 & mirna.cut == 3), which(methylation.cut == 1 & mirna.cut == 3), which(methylation.cut == 3 & mirna.cut == 1))
  groups = c(rep(1, sum(methylation.cut == 1 & mirna.cut == 1)), rep(2, sum(methylation.cut == 3 & mirna.cut == 3)), rep(3, sum(methylation.cut == 1 & mirna.cut == 3)), rep(4, sum(methylation.cut == 3 & mirna.cut == 1)))
  
  gd = data.frame(d$isoform.expression[quad], groups, stringsAsFactors = F)
  anova_summary = summary(aov(gd$d.isoform.expression.quad. ~ gd$groups))
  anova[i] = anova_summary[[1]][[5]][1]
}

map$kp_pval = kp_pval
map$cox_pval = cox_pval
map$t_test = t_test
map$anova = anova
write.table(map, "xena_results_iqr3_survival.csv", quote = T, col.names = T, row.names = F, sep=',')

