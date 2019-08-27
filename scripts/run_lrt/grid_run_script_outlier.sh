#!/bin/bash
#$ -N iqr3
#$ -cwd
#$ -l h_vmem=8G
#$ -l h_rt=3:00:00
#$ -t 1-400

Rscript run_lrt.R output_outlier_iqr1.5 IQR 1.5 $SGE_TASK_ID 400
