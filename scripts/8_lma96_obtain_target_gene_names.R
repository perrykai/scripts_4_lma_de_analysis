#' **Script:** `9_lma96_obtain_target_gene_names.R`
#' 
#' **Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/4_lma_de_analysis/scripts`
#' 
#' **Date:**  12/6/16
#' 
#' **Input File Directory:**  
#' 
#' 1. `/mnt/research/pigeqtl/analyses/microRNA/4_lma_de_analysis/`
#' 2. `/mnt/research/pigeqtl/analyses/microRNA/reference_sequences/targetscan_databases`
#' 
#' **Input File(s):** ``
#' 
#' 1. `9_targetscan_analysis_results.txt`
#' 2. `Gene_info.txt`
#' 
#' **Output File Directory:** 
#' 
#' `/mnt/research/pigeqtl/analyses/microRNA/4_lma_de_analysis/`
#' 
#' **Output File(s):** 
#' 
#' 1. `10_lma96_target_gene_transcript_ids_gene_symbols.txt`
#' 2. `11_lma96_target_transcript_ids.txt`
#' 
#' **Table of contents:**
#'
#' 1. [Objectives](#objectives)
#' 2. [Install libraries](#install-libraries)
#' 3. [Load data](#load-data)
#' 4. [Analysis](#analysis)
#' 5. [Visualize](#visualize)
#' 6. [Save data](#save-data)
#' 
#' ## Objectives
#' 
#' The objective of this script is to extract the transcript IDs from the unique targetscan results, and match their appropriate gene symbols for use in the DAVID/KEGG GO analysis
#' 
#' ## Install libraries
setwd("/mnt/research/pigeqtl/analyses/microRNA/4_lma_de_analysis/scripts/")

#' ## Load data
#' 
#' First, the target gene information
targetgenes<-read.table("../9_targetscan_analysis_results.txt", header=T, sep="\t", colClasses=c(rep("character", 2), rep("numeric", 6), rep("character", 3), rep("numeric", 2)))
dim(targetgenes)
head(targetgenes)

#' Next, the human gene info dataset (contains Gene ID, Gene Symbol, Description, etc)
geneinfo<-read.table("../../reference_sequences/targetscan_databases/Gene_info.txt", header=T, encoding="latin1",quote="",sep="\t", colClasses=c(rep("character",4),rep("numeric", 3)))
dim(geneinfo)
head(geneinfo)

#' ## Analysis
#' 
#' Eliminate the "6mer" site types, as these are the weakest likely targets:
targetgenes<-targetgenes[targetgenes$Site_type=="8mer-1a",]
#' First, I want to extract the transcript IDs from the gene targets (column 1)
transcriptids<-unique(targetgenes[,1])
#' How many unique trancript IDs are in the targets?
length(transcriptids)
head(transcriptids)

#' Use this vector as an index to extract the correct gene name from the gene_info dataset:
genesymbols<-geneinfo[match(transcriptids, geneinfo$Transcript.ID),c("Transcript.ID","Gene.symbol")]
dim(genesymbols)
head(genesymbols)

if(sum(genesymbols$Transcript.ID != transcriptids) != 0) stop ("Transcript ID match index did not work correctly")

#' ## Save data
write.table(genesymbols, file="../10_lma96_target_gene_transcript_ids_gene_symbols.txt", quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t")
write.table(genesymbols$Transcript.ID, file="../11_lma96_target_transcript_ids.txt", quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t")