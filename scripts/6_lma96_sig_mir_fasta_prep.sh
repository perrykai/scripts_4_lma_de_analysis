#' **Script:** `6_lma96_sig_mir_fasta_prep.sh`
#' 
#' **Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/4_lma_de_analysis/scripts`
#' 
#' **Date:**  12/01/16
#' 
#' **Input File Directory:**  
#' 
#' 1. `/mnt/research/pigeqtl/analyses/microRNA/4_lma_de_analysis`
#' 2. `/mnt/research/pigeqtl/analyses/microRNA/reference_sequences`
#' 
#' **Input File(s):** 
#' 
#' 1. `6_sig_mirna_names.txt`
#' 2. `ssc_mature_mir.fa`
#' 
#' **Output File Directory:** `/mnt/research/pigeqtl/analyses/microRNA/4_lma_de_analysis`
#' 
#' **Output File(s):** 
#' 
#' 1. `7_lma96_sig_mirna.txt`
#' 
#' **Table of contents:**
#'
#' 1. [Objectives](#objectives)
#' 2. [Save data](#save-data)
#' 
#' ## Objectives
#' 
#' The objective of this script is to extract the sequences of the significant DE mature miRNA from the reference file to compare seed sequences between pig and human miRNAs. 
#' This script will be executed in the script "5_targetscan_data_preparation.R"
#' 
#' ## Analysis
cd /mnt/research/pigeqtl/analyses/microRNA/4_lma_de_analysis
grep -A1 -f 6_sig_mirna_names.txt ../reference_sequences/ssc_mature_mir.fa | sed '/--/d' > 7_lma96_sig_mirna.fa

awk '{printf "%s%s",$0,(NR%2?FS:RS)}' 7_lma96_sig_mirna.fa > 7_lma96_sig_mirna.txt

rm 7_lma96_sig_mirna.fa