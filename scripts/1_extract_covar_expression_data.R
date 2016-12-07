#' **Script:** `1_extract_covar_expression_data.R`
#' 
#' **Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/4_lma_de_analysis/scripts`
#' 
#' **Date:**  11/8/16
#' 
#' **Input File Directory:**  
#' 
#' 1. `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/3_build_dge_object_for_eqtl/`
#' 2. `/mnt/research/pigeqtl/analyses/microRNA/1_preprocess_fastq_files/10_mirdeep2_core_quantify_predict_output/`
#' 3. `/mnt/research/pigeqtl/analyses/microRNA/1_preprocess_fastq_files/8_collapsed_fasta_output_expression_matrix/`
#' 
#' **Input File(s):** 
#' 
#' 1. `4_msuprp_mirna_gpdata_pheno_counts.Rdata`
#' 2. `miRNAs_expressed_all_samples_21_01_2016_t_20_01_10.csv`
#' 3. `config_for_mapper.txt`
#' 
#' **Output File Directory:** `/mnt/research/pigeqtl/analyses/microRNA/4_lma_de_analysis/`
#' 
#' **Output File(s):** `1_lma96_rounded_mean_mature_mirna_expression.Rdata`
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
#' The objective of this script is to extract the miRNA expression data pertaining to the 96 animals that will be used in the miRNA DE analysis for LMA extremes.
#' 
#' ## Install libraries
rm(list=ls())
library(plyr)

#' ## Load data
setwd("/mnt/research/pigeqtl/analyses/microRNA/4_lma_de_analysis/scripts")

load("../../2_mirna_characterization_expression/3_build_dge_object_for_eqtl/4_msuprp_mirna_gpdata_pheno_counts.Rdata")

#' Read in the csv of the miRNA expression profiles, using check.names default to correct column names with non-ASCII characters
rc<-read.csv("../../1_preprocess_fastq_files/10_mirdeep2_core_quantify_predict_output/miRNAs_expressed_all_samples_21_01_2016_t_20_01_10.csv", sep = "\t", header = TRUE, row.names=NULL)

#' Read in the config file, maintaining the characters in the 3-digit code names
configfile<-read.table("../../1_preprocess_fastq_files/8_collapsed_fasta_output_expression_matrix/config_for_mapper.txt", header = FALSE, sep = " ", row.names=NULL, colClasses = c('character','character'))
head(configfile)
dim(configfile)

#' Remove the "_express.fa" from each file name to leave the pig id:
configfile$V1<-gsub("_express.fa", "", configfile$V1)

#' Make filenames more informative:
colnames(configfile)<-c("pigid","code")
head(configfile)
colnames(configfile)

ls()

#' ## Analysis
#' 
#' First, extract the IDs from the animals meeting the selection criteria of "lma", loin muscle area.
names(final_MSUPRP_miRNA)
head(final_MSUPRP_miRNA$covar)

table(final_MSUPRP_miRNA$covar$selcrit)
final_MSUPRP_miRNA$covar$selcrit

lmapigs<-final_MSUPRP_miRNA$covar[final_MSUPRP_miRNA$covar$selcrit == "lma",]
dim(lmapigs)

table(lmapigs$Status)
table(lmapigs$sex)

ids<-lmapigs$id

head(rc[,1:5])

#' Set the name of the first column to "miRNA":
colnames(rc)[[1]]<-"miRNA"

#' Remove the "X" character from the beginning of each column:
colnames(rc)<-gsub("X", "", colnames(rc))
colnames(rc)[1:30]

#' ### 1. Extract the columns of mature read counts for each miRNA for each animal
lmaquant<-rc[,c(1,5:178)]
colnames(lmaquant)

dim(lmaquant)

head(lmaquant[1:8])


#' Take a subset of this data.frame for testing:
test<-lmaquant[1:20,1:8]

test[1:10,]

#' ### 2. Use the 'by' function to apply the function colMeans to the entire data frame of read counts:
#' (What this will do is go down the columns looking at the index of grouped miRNA names and take the average of the read counts for that miRNA)
#' The result of this will be a list containing the average read counts for each miRNA for each animal.
#' 
#' Example: by(data, index, function)
bytst<-by(test[,2:ncol(test)], test[,1], colMeans)
bytst[1:25]
#' Notice here that the rest of the miRNA names remain since the miRNA name is a factor, but since there is no data for them they are filled with NULL.
#' 

#' Apply the by function to the full dataframe:
meanrc<-by(lmaquant[,2:ncol(lmaquant)], lmaquant[,1], colMeans)

#' This should be 411 (the number of mature pig miRNAs in miRBase), meaning we have one expression profile for each mature miRNA:
length(meanrc)

meanrc[1:3]

#' 
#' ### 3. Transform the list output back into a data.frame using the plyr package to prepare for gblup function:

#' Example: ldply(.data, .fun, .id)
#' 
#' id = name of the index column (used if data is a named list). Pass NULL to avoid creation
#'      of the index column. For compatibility, omit this argument or pass "NA" to avoid converting the index column
#'      to a factor; in this case, ".id" is used as column name.

lma.dfmeanrc<-ldply(meanrc, fun=NULL, id=names(meanrc))

head(lma.dfmeanrc[1:8])

dim(lma.dfmeanrc)

#' These dimensions are what would be expected, because there are 411 mature sus scrofa miRNA sequences in miRBase,
#' and there are 174 animals in the analysis, plus the miRNA column.
#' 
#' Check that the correct miRNA name went with the correct data:
if (sum(names(meanrc)!=lma.dfmeanrc[,1]) != 0) stop ("miRNA names are not the same")

colnames(lma.dfmeanrc)[[1]]<-"miRNA"

if (sum(colnames(lma.dfmeanrc)!=colnames(lmaquant)) != 0) stop ("animal order not the same")

head(lma.dfmeanrc[,1:10])

#' 
#' Set first column of lma.dfmeanrc (miRNA ids) as the row.names:
rownames(lma.dfmeanrc)<-lma.dfmeanrc$miRNA

#' Eliminate column of row names:
lma.dfmeanrc<-lma.dfmeanrc[,-c(1)]
head(lma.dfmeanrc[,1:10])
dim(lma.dfmeanrc)

#' ### 4. Restore the pig IDs in place of the 3-digit codes as the column names of the data frame, for use with the gblup function of gwaR

head(configfile)

#' Now I need to substitute the 3-digit code with the pig IDs, ensuring the names stay in the correct order:
#' 
#' Use match function to find positional index and match column names:
#' 
#' The object lma.dfmeanrc has column names that need to be re-named. I have the config file which contains
#' the current column names and the desired column names. What I am doing in this code is re-ordering the config file 
#' based on where the config file "code" column matches the position of the lma.dfmeanrc object's column names, then having it return the corresponding value in column "pigid". 
#' 
#' So, when using match, need to have the first argument be the matrix/dataframe you want to change or match, and the second argument be what you want to index it by or match it against. 
#' 
#' "Where does [vector] match in [matrix]?" or "Match the column names of lma.dfmeanrc to the configfile "code" column, then return the corresponding pigid."
configfile[match(colnames(lma.dfmeanrc),configfile$code),"pigid"]

#' Assign the column names using match:
colnames(lma.dfmeanrc)<- configfile[match(colnames(lma.dfmeanrc),configfile$code),"pigid"]
head(lma.dfmeanrc[1:10])
dim(lma.dfmeanrc)

if (sum(colnames(lma.dfmeanrc)!=(configfile$pigid))!=0) stop ("match function did not work correctly")

#' ### 5. Round the mean read counts to the nearest integer for use with the voom function
lma.dfmeanrc[1:10,1:6]

lma.dfmeanrcround<-round(lma.dfmeanrc)

lma.dfmeanrcround[1:10,1:6]

#' The final matrix of rounded, mean read counts needs to have miRNA rownames and Animal ID colnames
head(rownames(lma.dfmeanrcround))

head(colnames(lma.dfmeanrcround))

if (sum(colnames(lma.dfmeanrcround)!= configfile$pigid)!= 0) stop ("rownames do not match pigid")

#' ### 6. Extract the 96 libraries from the 174 Bioo Scientific dataset matching the 96 animals with the "lma" selection criteria.
lma96.dfmeanrcround<-lma.dfmeanrcround[,ids]
dim(lma96.dfmeanrcround)
head(lma96.dfmeanrcround)

colnames(lma96.dfmeanrcround)
sum(colnames(lma96.dfmeanrcround) != ids)

#' Check that the data is the same between the subset and the lma.dfmeanrcround dataset:
for (i in colnames(lma96.dfmeanrcround)){
        print(all.equal(lma96.dfmeanrcround[,i], lma.dfmeanrcround[,i]))
}

#' How many miRNAs have 0 expression in this dataset?
sum(rowSums(lma96.dfmeanrcround)==0)

#' ## Visualize

#' ## Save data
#' What I am saving here is the rounded, average read counts, and the corresponding covariate data in an .Rdata object
save(lma96.dfmeanrcround, lmapigs, file="../1_lma96_covar_rounded_mean_mature_mirna_expression.Rdata")