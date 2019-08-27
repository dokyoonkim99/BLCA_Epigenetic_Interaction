library(stringr)
library(ggplot2)
library(interactions)

isoform_data = data.frame(fread("../data/qc_centered_data_xena/isoform_log_centered_data"), stringsAsFactors = F)
methylation_data = data.frame(fread("../data/qc_centered_data_xena/methylation_centered_data"), stringsAsFactors = F)
mirnaseq_data = data.frame(fread("../data/qc_centered_data_xena/mirna_log_centered_data"), stringsAsFactors = F)
clinical = data.frame(fread("../data/qc_centered_data_xena/clinical_data"), stringsAsFactors = F)
clinical$OS.time = clinical$OS.time /  30.417
clinical$surv = Surv(clinical$OS.time, clinical$OS)

map = read.csv(file = 'xena_results_iqr3_survival.csv', stringsAsFactors = F)
map = map[map$cox_pval < 0.05,]
isoform.name = "TGFBR3_ENST00000212355"
methylation.name = "cg08648138"
mirna.name = "hsa-let-7c-5p"

for (i in 1:nrow(map)) {
  m = as.vector(unlist(map[i,c(1:3)]))
  isoform.name = m[1]
  methylation.name = m[2]
  mirna.name = m[3]
  ret = paste(isoform.name, methylation.name, mirna.name, sep='_')
  
  isoform_row = which(isoform_data$V1 == isoform.name)
  methylation_row = which(methylation_data$V1 == methylation.name)
  mirna_row = which(mirnaseq_data$V1 == mirna.name)
  
  isoform.expression = unlist(isoform_data[isoform_row,-1])
  methylation = unlist(methylation_data[methylation_row,-1])
  mirna = unlist(mirnaseq_data[mirna_row,-1])
  
  d = data.frame(isoform.expression, methylation, mirna, clinical$gender, clinical$age_at_initial_pathologic_diagnosis ,clinical$ajcc_pathologic_tumor_stage, clinical$histological_grade, clinical$surv)
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

  dir.create(paste('figures/interaction',  sep=''), showWarnings = F)
  png(paste('figures/interaction/', ret, '.png', sep=''), width=7, height=5, units="in",res=300)
  #p = interact_plot(full, pred = 'methylation', modx = 'mirna', y.label = "Isoform", x.label = "Methylation", legend.main = "miRNA", 
   #                 main.title = paste(isoform.name, methylation.name, mirna.name, sep='  '))
  p = interact_plot(full, pred = 'mirna', modx = 'methylation', y.label = "Isoform", x.label = "miRNA", legend.main = "Methylation", 
                    main.title = paste(isoform.name, methylation.name, mirna.name, sep='  '), line.thickness = 5)
  #p = interact_plot(full, pred = 'mirna', modx = 'methylation')
  p = p + theme(plot.title = element_text(hjust = 0.5))
  plot(p)
  dev.off()
}
