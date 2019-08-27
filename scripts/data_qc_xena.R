library(data.table)
library(stringr)

# Isoform data 
isoform = fread(file = "../data/xena_data/tcga_RSEM_isoform_fpkm_BLCA")
isoform$sample = sapply(str_split(isoform$sample, "\\."), function(x){return(x[1])})
isoform_map = read.csv(file = "../data/xena_data/gencode.v23.annotation.transcript.probemap", sep =  '\t', stringsAsFactors = F)
isoform_map$id = sapply(str_split(isoform_map$id, "\\."), function(x){return(x[1])})
isoform_map$id_gene = paste( isoform_map$gene, isoform_map$id, sep = '_')
isoform_map = isoform_map[,c("id", "id_gene")]
isoform = merge(isoform_map, isoform, by.x = "id", by.y = "sample")
isoform = isoform[,-1]

colnames(isoform) = substr(colnames(isoform), 1, 15)
colnames(isoform) = gsub("-", ".", colnames(isoform))
rnames = isoform[,1]
isoform = isoform[,-1]
isoform = data.frame(isoform, stringsAsFactors = F)
isoform = as.data.frame(sapply(isoform, as.double))
rownames(isoform) = rnames

controls = grepl("....\\...\\.....\\.[1|2][0-9]*.", colnames(isoform))
colnames(isoform)[controls]
isoform = isoform[, !controls]
a = apply(isoform, 1, function(x){sum(x >= log2(0.1 + 0.001), na.rm = T)/length(x) > 0.5})
isoform = isoform[a,]

# Methylation data
methylation_file = "../data/xena_data/HumanMethylation450"
methylation = data.frame(fread(methylation_file), stringsAsFactors = F)
rnames = methylation[,1]
methylation = methylation[,-1]
methylation = as.data.frame(sapply(methylation, as.double))
row.names(methylation) = rnames
na = apply(methylation, 1, function(x){all(is.na(x))})
methylation = methylation[!na,]
controls = grepl("....\\...\\.....\\.[1|2][0-9]*.", colnames(methylation))
methylation = methylation[,!controls]

# miRNA
file = "../data/gdac.broadinstitute.org_BLCA.miRseq_Mature_Preprocess.Level_3.2016012800.0.0/BLCA.miRseq_mature_RPM_log2.txt"
mirnaseq = data.frame(fread(file))
colnames(mirnaseq) = substr(colnames(mirnaseq), 1, 15)
rnames = matrix(unlist(str_split(mirnaseq[,1], '\\|')), ncol=2, byrow=TRUE)[,1]
mirnaseq = mirnaseq[,-1]
mirnaseq = as.data.frame(sapply(mirnaseq, as.double))
rownames(mirnaseq) = rnames
na = apply(mirnaseq, 1, function(x){sum(is.na(x))/length(x) > 0.75})
mirnaseq = mirnaseq[!na,]
controls = grepl("....\\...\\.....\\.[1|2][0-9]*.", colnames(mirnaseq))
mirnaseq = mirnaseq[,!controls]


#clinical_pancan = read.csv(file = "Survival_SupplementalTable_S1_20171025_xena_sp", sep = '\t')
#clinical_pancan = clinical_pancan[clinical_pancan$cancer.type.abbreviation == "BLCA",]
#clinical_pancan = clinical_pancan[,c("sample", "OS", "OS.time")]

#clinical = clinical[, c("sampleID", "X_TIME_TO_EVENT", "X_EVENT")]
#merged = merge(clinical_pancan, clinical, by.x = "sample", by.y = "sampleID", all = T)
#diff = merged[merged$OS.time != merged$X_TIME_TO_EVENT, ]

# Clinical
clinical = read.csv(file = "../data/xena_data/Survival_SupplementalTable_S1_20171025_xena_sp", na.strings = c("na", "NA", "[Unknown]"), sep = '\t', stringsAsFactors = F)
clinical = clinical[clinical$cancer.type.abbreviation == "BLCA",]
# no DOB or discrepency
clinical = clinical[clinical$sample != "TCGA-HQ-A2OE-01", ]
# cancer stage NA
clinical = clinical[clinical$sample != "TCGA-FJ-A3Z9-01", ]
# no survival
clinical = clinical[clinical$sample != "TCGA.GV.A3QG.01", ]

clinical = clinical[, c('sample', 'age_at_initial_pathologic_diagnosis', 'OS.time', 'OS', 'gender', 
                        "ajcc_pathologic_tumor_stage", "histological_grade")]
clinical$sample = gsub("-", ".", clinical$sample)
clinical = clinical[!is.na(clinical$histological_grade),]
controls = grepl("....\\...\\.....\\.[1|2][0-9]*.", clinical$sample)
clinical = clinical[!controls,]

# table(is.na(clinical$pathologic_M) | is.na(clinical$pathologic_T) | is.na(clinical$pathologic_N))

# Common samples
common_samples = colnames(isoform)[colnames(isoform) %in% colnames(methylation) & colnames(isoform) %in% clinical$sample & colnames(isoform) %in% colnames(mirnaseq)]
isoform = isoform[,colnames(isoform) %in% common_samples]
methylation = methylation[,colnames(methylation) %in% common_samples]
mirnaseq = mirnaseq[,colnames(mirnaseq) %in% common_samples]
clinical = clinical[clinical$sample %in% common_samples,]

isoform = isoform[,order(colnames(isoform), common_samples)]
methylation = methylation[,order(colnames(methylation), common_samples)]
mirnaseq = mirnaseq[,order(colnames(mirnaseq), common_samples)]
clinical = clinical[order(clinical$sample, common_samples),]

mean = apply(isoform, 1, function(x){mean(x, na.rm =T)})
isoform = isoform - mean

mean = apply(methylation, 1, function(x){mean(as.numeric(x), na.rm =T)})
methylation = methylation - mean

mean = apply(mirnaseq, 1, function(x){mean(as.numeric(x), na.rm =T)})
mirnaseq = mirnaseq - mean

write.table(isoform, file = "../data/qc_centered_data_xena/isoform_log_centered_data", quote = F, sep = '\t')
write.table(methylation, file = "../data/qc_centered_data_xena/methylation_centered_data", quote = F, sep = '\t')
write.table(mirnaseq, file = "../data/qc_centered_data_xena/mirna_log_centered_data", quote = F, sep = '\t')
write.table(clinical, file = "../data/qc_centered_data_xena/clinical_data", quote = F, sep = '\t')

