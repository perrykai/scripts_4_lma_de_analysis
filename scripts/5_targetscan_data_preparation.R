#' **Script:** `5_targetscan_data_preparation.R`
#' 
#' **Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/4_lma_de_analysis/scripts`
#' 
#' **Date:**  11/29/16
#' 
#' **Input File Directory:**  
#' 
#' 1. `/mnt/research/pigeqtl/analyses/microRNA/reference_sequences/targetscan_databases`
#' 2. `/mnt/research/pigeqtl/analyses/microRNA/4_lma_de_analysis`
#' 
#' **Input File(s):** 
#' 
#' 1. `UTR_sequences_hsa.txt`
#' 2. `miR_Family_Info_hsa.txt`
#' 3. `5_lma96_voom_de_analysis_results.R`
#' 
#' **Output File Directory:** `/mnt/research/pigeqtl/analyses/microRNA/4_lma_de_analysis`
#' 
#' **Output File(s):** 
#' 
#' 1. `6_sig_mirna_names.txt`
#' 2. `7_lma96_sig_mirna.txt (output from script 6_lma96_sig_mir_fasta_prep.sh, executed in this script)`
#' 3. `8_lma96_final_sig_hsa_mirna.txt`
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
#' The objective of this script is to prepare the miRNA datasets for use in the targetscan target prediction analysis
#' 
#' The list of significant miRNAs identified in the DE analysis will be extracted and used as an index to subset the hsa reference datasets.
#'
#' This will allow the human-equivalent miRNAs to be queried against the human gene list for target prediction. 
#' 
#' Another point to consider would be to limit the UTR data to only those genes expressed in the 174 F2 pigs, using a gene list provided by Deborah.
#' 
#' ## Install libraries
#' ## Load data
setwd("/mnt/research/pigeqtl/analyses/microRNA/4_lma_de_analysis/scripts")
#' Load the results from the DE Analysis
load("../5_lma96_voom_de_analysis_results.Rdata")
#' Load the human miRNAs
hsamirna<-read.table("../../reference_sequences/targetscan_databases/miR_Family_Info_hsa.txt")
colnames(hsamirna)<-c("mirna_name", "seed_seq", "species_id")
head(hsamirna)
hsamirna$mirna_name<-as.character(hsamirna$mirna_name)

#' ## Analysis
#' 
#' Isolate the names of the miRNAs significant in the DE analysis
sum(toptab$adj.P.Val < 0.05)
sigmir<-toptab[toptab$adj.P.Val < 0.05,]
dim(sigmir)

pigmirnm<-as.character(sigmir$Name)
length(pigmirnm)
pigmir<-sub("ssc-", "", pigmirnm)
head(pigmir)
#' The following allows the significant pig miRNAs to be extracted correctly from the full dataset to compare seed sequences to human miR
pigmirnm2<-paste(pigmirnm, '$', sep="")
#' Write this to a .txt file right now, then execute the next script from this script using that file
write.table(pigmirnm2, file="../6_sig_mirna_names.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)
#' Change the permissions on the next script:
system("chmod 774 6_lma96_sig_mir_fasta_prep.sh")
#' Execute the next script so the analysis can continue
system("./6_lma96_sig_mir_fasta_prep.sh", wait=FALSE)

#' Need to remove the "3p", "5p" information to match between this dataset and the human dataset
#' (this is ugly, but it works - ask if there's a better way to do this, haha)

pigsub<-strsplit(pigmir, "-")
pigsub<-do.call(rbind, pigsub)[,1:2]
head(pigsub)
pigsub<-paste(pigsub[,1], pigsub[,2], sep="-")
head(pigsub)
length(pigsub)
#' How many unique miRNAs are in the dataset?
length(unique(pigsub))

#' Use these names to find matches against the human miRNA database
y<-list()
x<-list()
hsasub<-data.frame()
for(i in pigsub){
	y[[i]]<-grep(i, hsamirna$mirna_name)
	x[[i]]<-hsamirna[y[[i]],]
}
length(y)
length(x)

xdf<-do.call(rbind, x)
head(xdf)
dim(xdf)
#' How many miRNAs had a match in this step? (i.e., how many am I gonna have to go and find myself?)
length(unique(gsub("\\..*","",rownames(xdf))))

#' Instead of using the names to match the human miRNAs, use the seed sequences!
#' 
#' Load the known pig miRNA sequences corresponding to the significantly DE miRNAs
pigseq<-read.table("../7_lma96_sig_mirna.txt", col.names=c("name","seq"))
head(pigseq)
pigseq$name<-sub(">", "", pigseq$name)
#' Extract nt 2-8 (the miRNA seed sequence)
pigseq$seedseq<-substr(pigseq$seq, 2, 8)
head(pigseq)

#' Match the seed sequences between the pig seq and the human seq, returning the miRNA name of the matching human seq
pigseq$matchedseed<-xdf[match(pigseq$seedseq, xdf$seed_seq),"mirna_name"]
head(pigseq)
head(pigseq[,c(1,4)])
sum(is.na(pigseq$matchedseed))

#' This way matches the seed sequences between the full dataset of human miRNAs, not the filtered
pigseq$matchedseed2<-hsamirna[match(pigseq$seedseq, hsamirna$seed_seq), "mirna_name"]
head(pigseq[,c(1,5)])
sum(is.na(pigseq$matchedseed2))

sum(pigseq$matchedseed2 %in% pigseq$matchedseed)
dim(pigseq)
nrow(pigseq) - sum(pigseq$matchedseed2 %in% pigseq$matchedseed)
sum(is.na(pigseq$matchedseed)) - sum(is.na(pigseq$matchedseed2))
head(pigseq)

#' The next step is to use the name of the human miRNA sequences to subset the larger database for use with the targetscan algorithm
head(hsamirna)
finalhsamir<-hsamirna[match(pigseq$matchedseed2, hsamirna$mirna_name),]
finalhsamir<-finalhsamir[complete.cases(finalhsamir),]
dim(finalhsamir)
head(finalhsamir)

#' ## Save data
write.table(finalhsamir, file="../8_lma96_final_sig_hsa_mirna.txt", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")

