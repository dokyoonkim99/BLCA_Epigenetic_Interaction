import sys
import os

# generated by Seonggyun
combination_file = "../../data/07292019.miRNA-target-info.txt"
xena_methylation_gene_map = "../../data/xena_data/illuminaMethyl450_hg19_GPL16304_TCGAlegacy"

gene_cgid_map = {}
gene_hgnc_map = {}
with open(xena_methylation_gene_map) as xmf:
	next(xmf)
	for line in xmf:
		line = line.strip()
		parts = line.split('\t')
		genes = parts[1].split(',')
		for gene in genes:
			if gene not in gene_cgid_map:
				gene_cgid_map[gene] = set()
				gene_hgnc_map[gene] = gene
			gene_cgid_map[gene].add(parts[0])

# Accomodate old HGNC symbols used in TCGA data from xena
# gene_hgnc_map["SZRD1"] = "C1orf144"
# gene_hgnc_map["THEMIS2"] = "C1orf38"
# gene_hgnc_map["EXO5"] = "DEM1"
# gene_hgnc_map["TRMT13"] = "CCDC76"
# gene_hgnc_map["RTCA"] = "RTCD1"
# gene_hgnc_map["PIFO"] = "C1orf88"
gene_hgnc_map["01-Mar"] = "MARCH1"
# gene_hgnc_map["SPRTN"] = "C1orf124"
# gene_hgnc_map["COX20"] = "FAM36A"
# gene_hgnc_map["EMC1"] = "KIAA0090"
# gene_hgnc_map["STPG1"] = "C1orf201"
# gene_hgnc_map["AUNIP"] = "C1orf135"
# gene_hgnc_map["LAMTOR5"] = "HBXIP"
# gene_hgnc_map["MPC2"] = "BRP44"
# gene_hgnc_map["CCSAP"] = "C1orf96"
# gene_hgnc_map["FAM228A"] = "C2orf84"
# gene_hgnc_map["ATRAID"] = "C2orf28"
# gene_hgnc_map["NDUFAF7"] = "C2orf56"
# gene_hgnc_map["NABP1"] = "OBFC2A"
gene_hgnc_map["02-Sep"] = "SEPT2"
# gene_hgnc_map["TRABD2A"] = "C2orf89"
gene_hgnc_map["10-Sep"] = "SEPT10"
# gene_hgnc_map["KANSL1L"] = "C2orf67"
# gene_hgnc_map["CMSS1"] = "C3orf26"
# gene_hgnc_map["TRMT10C"] = "RG9MTD1"
# gene_hgnc_map["NXPE3"] = "FAM55C"
# gene_hgnc_map["COPG1"] = "COPG"
# gene_hgnc_map["EMC3"] = "TMEM111"
# gene_hgnc_map["ELP6"] = "C3orf75"
# gene_hgnc_map["EOGT"] = "C3orf64"
# gene_hgnc_map["GCSAM"] = "GCET2"
# gene_hgnc_map["MSANTD1"] = "C4orf44"
# gene_hgnc_map["TENM3"] = "ODZ3"
# gene_hgnc_map["SETD9"] = "C5orf35"
# gene_hgnc_map["ARHGEF28"] = "RGNEF"
# gene_hgnc_map["SPDL1"] = "CCDC99"
gene_hgnc_map["03-Mar"] = "MARCH3"
# gene_hgnc_map["ABRACL"] = "C6orf115"
# gene_hgnc_map["GINM1"] = "C6orf72"
# gene_hgnc_map["CCDC170"] = "C6orf97"
# gene_hgnc_map["BLOC1S5"] = "MUTED"
# gene_hgnc_map["OARD1"] = "C6orf130"
# gene_hgnc_map["PTCHD4"] = "C6orf138"
# gene_hgnc_map["FAXC"] = "C6orf168"
# gene_hgnc_map["SLC18B1"] = "C6orf192"
# gene_hgnc_map["MPC1"] = "BRP44L"
# gene_hgnc_map["AP5Z1"] = "KIAA0415"
# gene_hgnc_map["FAM221A"] = "C7orf46"
gene_hgnc_map["07-Sep"] = "SEPT7"
# gene_hgnc_map["TMEM248"] = "C7orf42"
# gene_hgnc_map["LAMTOR4"] = "C7orf59"
# gene_hgnc_map["CPED1"] = "C7orf58"
# gene_hgnc_map["COA1"] = "C7orf44"
# gene_hgnc_map["MCMDC2"] = "C8orf45"
# gene_hgnc_map["NDUFAF6"] = "C8orf38"
# gene_hgnc_map["DCSTAMP"] = "TM7SF4"
# gene_hgnc_map["EMC2"] = "TTC35"
# gene_hgnc_map["THEM6"] = "C8orf55"
# gene_hgnc_map["CCDC171"] = "C9orf93"
# gene_hgnc_map["FOCAD"] = "KIAA1797"
# gene_hgnc_map["TRMT10B"] = "RG9MTD3"
# gene_hgnc_map["ERCC6L2"] = "C9orf102"
# gene_hgnc_map["MSANTD3"] = "C9orf30"
# gene_hgnc_map["NTMT1"] = "METTL11A"
# gene_hgnc_map["CACFD1"] = "C9orf7"
# gene_hgnc_map["RABL6"] = "C9orf86"
# gene_hgnc_map["PLGRKT"] = "C9orf46"
# gene_hgnc_map["ARHGEF39"] = "C9orf100"
# gene_hgnc_map["FAM221B"] = "C9orf128"
# gene_hgnc_map["SLC25A51"] = "MCART1"
# gene_hgnc_map["NMRK1"] = "C9orf95"
# gene_hgnc_map["TMEM254"] = "C10orf57"
gene_hgnc_map["05-Mar"] = "MARCH5"
# gene_hgnc_map["R3HCC1L"] = "C10orf28"
# gene_hgnc_map["PLEKHS1"] = "C10orf81"
# gene_hgnc_map["SKIDA1"] = "C10orf140"
gene_hgnc_map["08-Mar"] = "MARCH8"
# gene_hgnc_map["MSS51"] = "ZMYND17"
# gene_hgnc_map["ARL14EP"] = "C11orf46"
# gene_hgnc_map["KIAA1549L"] = "C11orf41"
# gene_hgnc_map["VPS51"] = "C11orf2"
# gene_hgnc_map["AAMDC"] = "C11orf67"
# gene_hgnc_map["TENM4"] = "ODZ4"
# gene_hgnc_map["MSANTD4"] = "KIAA1826"
# gene_hgnc_map["KIAA1551"] = "C12orf35"
# gene_hgnc_map["ASIC1"] = "ACCN2"
# gene_hgnc_map["METTL25"] = "C12orf26"
# gene_hgnc_map["HECTD4"] = "C12orf51"
# gene_hgnc_map["VWA8"] = "KIAA0564"
# gene_hgnc_map["TEX30"] = "C13orf27"
# gene_hgnc_map["AP5M1"] = "MUDENG"
# gene_hgnc_map["PCNXL4"] = "C14orf135"
# gene_hgnc_map["CCDC176"] = "C14orf45"
# gene_hgnc_map["ZC2HC1C"] = "FAM164C"
# gene_hgnc_map["GSKIP"] = "C14orf129"
# gene_hgnc_map["APOPT1"] = "C14orf153"
# gene_hgnc_map["DTD2"] = "C14orf126"
# gene_hgnc_map["L3HYPDH"] = "C14orf149"
# gene_hgnc_map["VIPAS39"] = "C14orf133"
# gene_hgnc_map["NRDE2"] = "C14orf102"
# gene_hgnc_map["TICRR"] = "C15orf42"
# gene_hgnc_map["EMC7"] = "C15orf24"
# gene_hgnc_map["KATNBL1"] = "C15orf29"
# gene_hgnc_map["PKM"] = "PKM2"
# gene_hgnc_map["KDM8"] = "JMJD5"
# gene_hgnc_map["CNEP1R1"] = "TMEM188"
# gene_hgnc_map["USB1"] = "C16orf57"
# gene_hgnc_map["MSRB1"] = "SEPX1"
# gene_hgnc_map["CMC2"] = "C16orf61"
# gene_hgnc_map["EMC8"] = "COX4NB"
# gene_hgnc_map["EFCAB13"] = "C17orf57"
gene_hgnc_map["09-Sep"] = "SEPT9"
# gene_hgnc_map["KANSL1"] = "KIAA1267"
gene_hgnc_map["04-Sep"] = "SEPT4"
# gene_hgnc_map["HID1"] = "C17orf28"
# gene_hgnc_map["LDLRAD4"] = "C18orf1"
gene_hgnc_map["02-Mar"] = "MARCH2" 
# gene_hgnc_map["SMIM7"] = "C19orf42"
# gene_hgnc_map["DNAAF3"] = "C19orf51"
# gene_hgnc_map["AAR2"] = "C20orf4"
# gene_hgnc_map["RTFDC1"] = "C20orf43"
# gene_hgnc_map["GID8"] = "C20orf11"
# gene_hgnc_map["SLC52A3"] = "C20orf54"
# gene_hgnc_map["HELZ2"] = "PRIC285"
gene_hgnc_map["05-Sep"] = "SEPT5"
gene_hgnc_map["03-Sep"] = "SEPT3"
# gene_hgnc_map["DESI1"] = "PPPDE2"
# gene_hgnc_map["DENND6B"] = "FAM116B"
# gene_hgnc_map["EPPIN"] = "SPINLW1"
# gene_hgnc_map["NUGGC"] = "C8orf80"
# gene_hgnc_map["TMEM246"] = "C9orf125"

ignore_set = set()
# ignore_set = set(["CFHR3", "GPR75-ASB3", "SERF1A", "IL33", "KLRC3", "GAS6", "CHP1", "01-Sep", "CCL3L1", "FAM227A", "LRRC19"])

combinations_set = set() 
with open(combination_file) as cf:
	for line in cf:
		line = line.strip()
		parts = line.split('\t')
		isoforms = parts[3].replace('"', '')
		isoforms = isoforms.split(',')
		for isoform in isoforms:
			if parts[0] not in ignore_set:
				for cgid in gene_cgid_map[gene_hgnc_map[parts[0]]]:
					combinations_set.add(parts[0] + '_' + isoform + '\t' + cgid + '\t' + parts[1] + '\t' + parts[0])


for line in combinations_set:
	print line
