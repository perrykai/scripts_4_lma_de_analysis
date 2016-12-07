#!/bin/sh  -login
#PBS -l nodes=1:ppn=1,walltime=03:00:00,mem=5Gb
#PBS -N targetscan_lma96_de_hsa_query
#PBS -j oe
#PBS -o /mnt/research/pigeqtl/analyses/microRNA/OutputsErrors/
#PBS -m abe
#PBS -M perrykai@msu.edu

#' **Script:** `7_run_targetscan_analysis.sh`
#' 
#' **Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/4_lma_de_analysis/scripts`
#' 
#' **Date:**  12/5/16
#' 
#' **Input File Directory:**  
#' 
#' 1. `/mnt/research/pigeqtl/analyses/microRNA/4_lma_de_analysis/`
#' 2. `/mnt/research/pigeqtl/analyses/microRNA/reference_sequences/targetscan_databases/`
#' 
#' **Input File(s):** 
#' 
#' 1. `8_lma96_final_sig_hsa_mirna.txt`
#' 2. `UTR_sequences_hsa.txt`
#' 
#' **Output File Directory:** 
#' 
#' 1. `/mnt/research/pigeqtl/analyses/microRNA/4_lma_de_analysis/`
#' 
#' **Output File(s):** 
#' 
#' 1. `9_targetscan_analysis_results.txt`
#' 
#' **Table of contents:**
#'
#' 1. [Objectives](#objectives)
#' 
#' ## Objectives
#' 
#' The objective of this script is to run the targetscan analysis using the 128 hsa miRNAs that had matching seed sequences to the significant DE miRNAs.
#' 

cd /mnt/research/pigeqtl/analyses/microRNA/4_lma_de_analysis/scripts

/mnt/home/perrykai/targetscan_70/targetscan_70.pl /mnt/research/pigeqtl/analyses/microRNA/4_lma_de_analysis/8_lma96_final_sig_hsa_mirna.txt /mnt/research/pigeqtl/analyses/microRNA/reference_sequences/targetscan_databases/UTR_sequences_hsa.txt /mnt/research/pigeqtl/analyses/microRNA/4_lma_de_analysis/9_targetscan_analysis_results.txt

qstat -f ${PBS_JOBID}