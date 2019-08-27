library(survival)
library(stringr)
library(ggplot2)
library(GGally)
library(data.table)
source('multiplot.R')

outlier_type = "IQR" 
outlier_threshold = 3

isoform_data = data.frame(fread("../data/qc_centered_data_xena/isoform_log_centered_data"), stringsAsFactors = F)
methylation_data = data.frame(fread("../data/qc_centered_data_xena/methylation_centered_data"), stringsAsFactors = F)
mirnaseq_data = data.frame(fread("../data/qc_centered_data_xena/mirna_log_centered_data"), stringsAsFactors = F)
clinical = data.frame(fread("../data/qc_centered_data_xena/clinical_data"), stringsAsFactors = F)
clinical$OS.time = clinical$OS.time /  30.417
clinical$surv = Surv(clinical$OS.time, clinical$OS)

map = read.csv(file = 'xena_results_iqr3_survival.csv', stringsAsFactors = F)
map = map[map$cox_pval < 0.05,]
sum(map$cox_pval < 0.05)
results = vector(length = nrow(map))
count = 0

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

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

plot_combined = function(data, num) {
  d = data
  methylation.breaks = c(quantile(d$methylation, probs = seq(0, 1, by = 0.333333333333333333333333333333333333333)))
  methylation.cut = cut(d$methylation, breaks=methylation.breaks, labels=c('1','2','3'), include.lowest=TRUE)
  
  mirna.breaks = c(quantile(d$mirna, probs = seq(0, 1, by = 0.333333333333333333333333333333333333333)))
  mirna.cut = cut(d$mirna, breaks=mirna.breaks, labels=c('1','2','3'), include.lowest=TRUE)
  
  quad = c(which(methylation.cut == 1 & mirna.cut == 1), which(methylation.cut == 3 & mirna.cut == 3))
  groups = c(rep(1, sum(methylation.cut == 1 & mirna.cut == 1)), rep(2, sum(methylation.cut == 3 & mirna.cut == 3)))
  
  d$V3 = NA
  d[which(methylation.cut == 1 & mirna.cut == 1), ncol(d)] = 1
  d[which(methylation.cut == 3 & mirna.cut == 3), ncol(d)] = 2
  
  surv_d = data[quad,]
  surv_d$groups = groups
  fit  = survfit(clinical.surv ~ groups, data = surv_d)
  sdiff = survdiff(clinical.surv ~ groups, data = surv_d)
  #p.val = 1 - pchisq(sdiff$chisq, length(sdiff$n) - 1)
  
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
  p.val = cox_summary$coefficients[1,5]
  
  d$V3[is.na(d$V3)] = 3
  
  gg_color_custom = gg_color_hue(4)
  gg_color_custom = gg_color_custom[c(3,1,2,4)]
  qp = ggplot(data = d, aes(methylation, mirna )) + geom_point(aes(colour = factor(V3)))
  qp = qp + geom_hline(yintercept = mirna.breaks[2])
  qp = qp + geom_hline(yintercept = mirna.breaks[3])
  qp = qp + geom_vline(xintercept = methylation.breaks[2])
  qp = qp + geom_vline(xintercept = methylation.breaks[3])
  #qp = qp + theme(legend.justification = c(0, 1), legend.position=c(0.01, 0.99))
  qp = qp + guides(linetype = FALSE) + scale_colour_manual(name = NULL, breaks = c(1,2,3), values = c(gg_color_custom[c(1,2)], "darkgrey"), labels = c(paste('Group 1 (L/L): ', fit$n[1], sep =''), paste('Group 2 (H/H): ', fit$n[2], sep =''), paste("remaining:",sum(d$V3 == 3))))
  qp = qp + ggtitle(paste('Isoform:', isoform.name, '  Methylation:', methylation.name, '\nmiRNA:', mirna.name, '\n', sep=""))
  qp = qp + xlab('Methylation\n\na) miRNA vs methylation showing distribution in quantiles 3x3.') + ylab('miRNA') + theme(plot.title = element_text(hjust = 0.5))
  
  
  s = ggsurv(fit, order.legend = F) + guides(linetype = FALSE)+ scale_colour_manual(name = NULL, values = gg_color_custom[c(1,2)], breaks = c(1,2),labels = c(paste('Group 1 (L/L): ', fit$n[1], sep =''), paste('Group 2 (H/H): ', fit$n[2], sep ='')))
  s = s + theme(legend.justification = c(0,0), legend.position=c(0.7, 0.77))
  s = s + annotate("text", x= 0.85 * max(fit$time), y=0.75, label= paste(paste('p-value: ', format(p.val, scientific=T, digits = 3), sep='')))
  s = s + xlab(paste('Time\n\nb) Survival analysis between LL and HH group.\n', sep=""))
  s = s + theme(plot.margin = unit(c(0,0.2,0,0), "cm"))
  
  group_names = c("Group LL", "Group HH", "Group LH", "Group HL")
  gd = data.frame(data$isoform.expression[quad],  group_names[groups])
  colnames(gd) = c("isoform_expression", "group")
  ttest.pval = t.test(gd$isoform_expression ~ gd$group)
  p = ggplot(gd, aes(x = factor(group, levels = group_names), y = isoform_expression, fill = factor(group, levels = group_names)))
  p = p + geom_boxplot(show.legend = F) + scale_fill_manual(values = gg_color_custom[c(1,2)])
  p = p + annotate("text", x = 1.5, y = max(gd$isoform_expression), label= paste(paste('p-value: ', format(ttest.pval$p.value, scientific=T, digits=3), sep='')))
  p = p + labs(y = "Isoform Expressoin", x = paste("\nc) Isoform expression boxplot for ", isoform.name ," isoform showing \nhigh expression for Group 1 over Group 2 samples.\n", sep = ""))
  p = p + theme(plot.margin = unit(c(0.6,0.2,0,0), "cm"))
  
  quad = c(which(methylation.cut == 1 & mirna.cut == 1), which(methylation.cut == 3 & mirna.cut == 3), which(methylation.cut == 1 & mirna.cut == 3), which(methylation.cut == 3 & mirna.cut == 1))
  groups = c(rep(1, sum(methylation.cut == 1 & mirna.cut == 1)), rep(2, sum(methylation.cut == 3 & mirna.cut == 3)), rep(3, sum(methylation.cut == 1 & mirna.cut == 3)), rep(4, sum(methylation.cut == 3 & mirna.cut == 1)))
  
  
  #dir.create(paste('figures/isoform_exp',  sep=''), showWarnings = F)
  gd = data.frame(data$isoform.expression[quad], group_names[groups], stringsAsFactors = F)
  #gd$group = factor(gd$group, levels = group_names)
  colnames(gd) = c("isoform_expression", "group")
  ttest.pval = summary(aov(gd$isoform_expression ~ gd$group))
  ttest.pval = ttest.pval[[1]][[5]][1]
  ag = ggplot(gd, aes(x = factor(group, levels = group_names), y = isoform_expression, fill = factor(group, levels = group_names)))
  ag = ag + geom_boxplot(show.legend = F) + scale_fill_manual(values = gg_color_custom)
  ag = ag + annotate("text", x = 2.5, y = 0.975 * max(gd$isoform_expression), label = paste("\nanova p-value :", format(ttest.pval, scientific=T, digits=3)))
  ag = ag + labs(y = "Isoform Expression", x = "\nd) Isoform expression boxplot for all groups")
  #ag = ag + ggtitle(paste('Isoform:', isoform.name, '   Methylation:', methylation.name, '   miRNA:', mirna.name, '\n', sep=""))
  ag = ag + theme(plot.title = element_text(hjust = 0.5)) 
  #plot(p)
  
  surv_d = data[quad,]
  surv_d$groups = group_names[groups]
  surv_d$groups = factor(surv_d$groups, levels = group_names)
  
  #surv_d = data.frame(surv[quad], groups)
  sdiff = survdiff(clinical.surv ~ groups, data = surv_d)
  fit  = survfit(clinical.surv ~ groups, data = surv_d)
  p.val = 1 - pchisq(sdiff$chisq, length(sdiff$n) - 1)
  
  as = ggsurv(fit, order.legend = F) + guides(linetype = FALSE)+ scale_colour_manual(name = group_names, values = gg_color_custom , breaks = group_names, labels = c(paste('Group 1 (L/L): ', fit$n[1], sep =''), paste('Group 2 (H/H): ', fit$n[2], sep = ''), paste('Group 3 (L/H): ', fit$n[3], sep = ''), paste('Group 4 (H/L): ', fit$n[4], sep ='')))
  #as = ggsurv(fit, order.legend = F, size.est = 3, plot.cens = F) + guides(linetype = FALSE)+ scale_colour_manual(name = group_names, values = gg_color_custom , breaks = group_names, labels = c(paste('Group 1 (L/L): ', fit$n[1], sep =''), paste('Group 2 (H/H): ', fit$n[2], sep = ''), paste('Group 3 (L/H): ', fit$n[3], sep = ''), paste('Group 4 (H/L): ', fit$n[4], sep ='')))
  as = as + annotate("text", x=max(fit$time)-(1/6)*max(fit$time), y=0.95, label= paste(paste('p-value: ', format(p.val, scientific=T, digits = 3), sep='')))
  as = as + xlab('Time(months)\n\ne) Survival analysis across 4 groups.')  #+ ggtitle(paste('Isoform:', isoform.name, '   Methylation:', methylation.name, '   miRNA:', mirna.name, '\n', sep=""))
  as = as + theme(legend.justification = c(0,0), legend.position=c(0, 0)) + theme(legend.title = element_blank())
  #as = as + theme(legend.position = "none")
  
  layout = matrix(c(0, 1, 1, 0 , 2, 2, 3, 3, 4, 4, 5, 5), nrow = 3, byrow = TRUE)
  
  png(paste('figures/combined/', 'combined_', num, '.png', sep=''), width=11.5, height=13, units="in",res=300)
  #png(paste('figures/combined/', 'combined_', num, '.png', sep=''), width=9.5, height=11, units="in",res=300)
  multiplot(qp, s, p, ag, as, layout=layout)
  dev.off()
}

isoform.name = "SGCD_ENST00000435422"
methylation.name = "cg19748027"
mirna.name = "hsa-miR-409-3p"

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
  
  data = data.frame(isoform.expression, methylation, mirna, clinical$gender, clinical$age_at_initial_pathologic_diagnosis ,clinical$ajcc_pathologic_tumor_stage, clinical$histological_grade, clinical$surv)
  data = data[!apply(data, 1, anyNA),]
  data = data[!is_outlier(data$isoform.expression, outlier_type, outlier_threshold) & !is_outlier(data$methylation, outlier_type, outlier_threshold) & !is_outlier(data$mirna, outlier_type, outlier_threshold), ]
  
  dir.create(paste('figures/combined',  sep=''), showWarnings = F)
  plot_combined(data, ret)
}

