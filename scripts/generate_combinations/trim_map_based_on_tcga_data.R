library(data.table)
map = read.csv(file = "transcript_methylation_mirna_combinations", sep = '\t', header = F, stringsAsFactors = F)

transcript = data.frame(fread("../../data/qc_centered_data/isoform_log_centered_data"), stringsAsFactors = F)
methylation = data.frame(fread("../../data/qc_centered_data/methylation_centered_data"), stringsAsFactors = F)
mirna = data.frame(fread("../../data/qc_centered_data/mirna_log_centered_data"), stringsAsFactors = F)
map = map[map$V1 %in% transcript$V1, ]
map = map[map$V2 %in% methylation$V1,]
map = map[map$V3 %in% mirna$V1,]

write.table(map, file = "transcript_methylation_mirna_combinations_tcga", sep = '\t', col.names = F, row.names = F, quote = F)
