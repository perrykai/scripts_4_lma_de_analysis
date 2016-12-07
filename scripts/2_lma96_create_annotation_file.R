#' **Script:** `2_lma96_create_annotation_file.R`
#' 
#' **Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/4_lma_de_analysis/scripts`
#' 
#' **Date:**  11/10/16
#' 
#' **Input File Directory:**  
#' 
#' 1. `/mnt/research/pigeqtl/analyses/microRNA/reference_sequences/`
#' 2. `/mnt/research/pigeqtl/analyses/microRNA/4_lma_de_analysis/`
#' 
#' **Input File(s):** 
#' 
#' 1. `scc.gff3`
#' 2. `1_lma96_rounded_mean_mature_mirna_expression.Rdata`
#' 
#' **Output File Directory:** `/mnt/research/pigeqtl/analyses/microRNA/4_lma_de_analysis/`
#' 
#' **Output File(s):** 
#'  
#' 1. `2_lma96_precursor_mirna_annotation.Rdata`
#' 2. `3_lma96_mature_mirna_annotation.Rdata`
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
#' The objective of this script is to build the annotation file for the lma96 miRNA dataset, for use in the LMA DE analysis.
#' 
#' ## Install libraries
#' 
library(rtracklayer)

#' ## Load data
rm(list=ls())

setwd("/mnt/research/pigeqtl/analyses/microRNA/4_lma_de_analysis/scripts/")

load("../1_lma96_covar_rounded_mean_mature_mirna_expression.Rdata")


GFFfile <- "../../reference_sequences/ssc.gff3"
gff <- import.gff(GFFfile, version = "3")

slotNames(gff)

gff

#' ## Analysis
#' 
#' Create a data frame of annotation data including the following columns:
#' 
#' Name, chr0, start, end, width, strand, type, Alias, Derives_from
#' 
#' Extract specific fields of information, build annotation data frame:
mirannot <- data.frame(unlist(elementMetadata(gff)$Name),
				  as.data.frame(seqnames(gff)),
				  as.data.frame(ranges(gff)),
                  as.data.frame(strand(gff)),
                  unlist(elementMetadata(gff)$type),
                  unlist(elementMetadata(gff)$Alias),
                  unlist(elementMetadata(gff)$Derives_from)
                  )


colnames(mirannot) <- c("Name", 
					  "chr0",
					  "start", "end", "width", 
					  "strand", 
					  "type", 
					  "Alias", 
					  "Derives_from"
					  )


#' Check the class of each column in the data frame:
colclass<-NULL

for (i in colnames(mirannot)) {
   colclass<-c(colclass,class(mirannot[,i]))
}

colclass

head(mirannot)

#' ### Isolate precursor miRNA annotation/assembly information
sum(mirannot$type=="miRNA_primary_transcript")

precursor.mirannot <- mirannot[(mirannot$type=="miRNA_primary_transcript"),]

dim(precursor.mirannot)
head(precursor.mirannot)

if (sum(precursor.mirannot$type == "miRNA") != 0) stop ("mature miRNA in primary miRNA annotation")

#' ### Isolate mature miRNA annotation/assembly information
sum(mirannot$type==("miRNA"))
mature.mirannot <- mirannot[(mirannot$type=="miRNA"),]

#' Order the rows by mature miRNA Name
mature.mirannot <- mature.mirannot[order(mature.mirannot$Name),]

dim(mature.mirannot)
head(mature.mirannot)
tail(mature.mirannot)

if (sum(mature.mirannot$type == "miRNA_primary_transcript") != 0) stop ("primary miRNA in mature miRNA annotation")

sum(table(mature.mirannot$Name)>1)

#' ### Extract the duplicated mature miRNAs, order by miRNA name:
dups<-data.frame(name=mature.mirannot$Name,dup1=duplicated(mature.mirannot$Name),dup2=duplicated(mature.mirannot$Name,fromLast=TRUE))
#' The "duplicated" function returns "TRUE" for the second and higher occurrence of a given element.
#' 
#' Executed "duplicated" on the data frame with normal order of rows, then reversed order of rows, then combined the logical index vectors.
head(dups)

rs<-rowSums(dups[,2:3])
#' As a check, calculated the rowSums of the data.frame and observed where the overlaps were.
#' 
#' Then, to calculate the total number of assembled duplicated miRNAs, subtract the overlaps from the total sum:
rs
sum(rs)
sum(rs>1)
sum(rs)-(sum(rs>1))


dup.maturemir <- mature.mirannot[(duplicated(mature.mirannot$Name) | duplicated(mature.mirannot$Name, fromLast=TRUE)),]
dim(dup.maturemir)
head(dup.maturemir)
tail(dup.maturemir)


#' ### Extract the precursor IDs for duplicated mature miRNAs
dup_derive <- as.data.frame(tapply(as.character(dup.maturemir$Derives_from), as.character(dup.maturemir$Name), paste))
dup_derive$Name <- rownames(dup_derive)
colnames(dup_derive) <- c("mult_precursors", "Name")

head(dup_derive)
dim(dup_derive)

if ((nrow(dup_derive) == length(unique(dup.maturemir$Name))) != "TRUE") stop ("extraction of duplicated mature miRNAs did not work correctly")

#' ### Build a data.frame containing the mature miRNAs, without duplicate information (all mature miRNAs present only once)
#'
sum(duplicated(mature.mirannot$Name)==FALSE)
single.mature.mirannot<-mature.mirannot[duplicated(mature.mirannot$Name)==FALSE,]

dim(single.mature.mirannot)
head(single.mature.mirannot)

if (sum(duplicated(single.mature.mirannot$Name))!=0) stop ("nonunique miRNA present in single.mature.mirannot")

#' ### Merge the single.mature.mirannot and the dup_derive data frames to obtain the multiple precursor sequences for a duplicated mature miRNA:
merged.single.mature.mirannot<-merge(single.mature.mirannot, dup_derive, by = "Name", all.x = TRUE)
rownames(merged.single.mature.mirannot)<-merged.single.mature.mirannot$Name

dim(merged.single.mature.mirannot)
head(merged.single.mature.mirannot)

if (sum(merged.single.mature.mirannot$Name != single.mature.mirannot$Name) != 0) stop ("merged single miRNA names not the same as single miRNA names")
if (sum(rownames(merged.single.mature.mirannot) != merged.single.mature.mirannot$Name) != 0) stop ("merged.single.mature.mirannot rownames not the same as its miRNA column")

#' Determine which miRNAs lack assembly information and extract them:
summary(match(rownames(lma96.dfmeanrcround), merged.single.mature.mirannot$Name))
length(rownames(lma96.dfmeanrcround))
length(merged.single.mature.mirannot$Name)

noannot<-rownames(lma96.dfmeanrcround)[is.na(match(rownames(lma96.dfmeanrcround), merged.single.mature.mirannot$Name))]

length(noannot)

#' The following 49 miRNAs lack assembly information:
noannot

noannot<-data.frame(Name=noannot, chr0=NA, start=NA, end=NA, width=NA, strand=NA, type=NA, Alias=NA, Derives_from=NA, mult_precursors=NA)
dim(noannot)
head(noannot)

#' Add the miRNAs lacking assembly information into the annotation data frame:
#' 
#' Extract the miRNA names from the count matrix:
count.mirnas<-as.data.frame(rownames(lma96.dfmeanrcround))
colnames(count.mirnas)<- "Name"
dim(count.mirnas)
head(count.mirnas)

#' Merge the annotation file with the miRNA names from the count matrix:
lma96.total.mature.annot<-merge(count.mirnas, merged.single.mature.mirannot, by = "Name", all.x = TRUE)
rownames(lma96.total.mature.annot)<-lma96.total.mature.annot$Name

dim(lma96.total.mature.annot)
head(lma96.total.mature.annot)

#' Check that the unassembled miRNAs bound correctly:
match(noannot$Name, lma96.total.mature.annot$Name)

lma96.total.mature.annot[match(noannot$Name, lma96.total.mature.annot$Name),]

if (sum(rownames(lma96.dfmeanrcround)!=rownames(lma96.total.mature.annot))!=0) stop ("rownames of count mx and annotation df not equal")

#' Combine the two columns of precursor information into one:
multp<-sapply(lma96.total.mature.annot[,"mult_precursors"],paste,collapse=",")
singlep<-as.character(lma96.total.mature.annot[,"Derives_from"])
multp[multp=="NA"]<-singlep[multp=="NA"]
head(multp)

lma96.total.mature.annot$Precursors<-as.character(multp)

lma96.total.mature.annot2<-lma96.total.mature.annot[,-c(9,10)]

if (sum(rownames(lma96.total.mature.annot2)!=count.mirnas$Name)!=0) stop ("rownames of count mx and annotation df not equal")
if (sum(rownames(lma96.total.mature.annot2)!=rownames(lma96.dfmeanrcround))!=0) stop ("rownames of count mx and annotation df not equal")

head(lma96.total.mature.annot2)

#' ## Save data
#' 
#' Save both annotation of mature miRNAs for the eQTL analysis, and annotation of precursor sequences for reference. 

save(lma96.total.mature.annot2, file = "../3_lma96_mature_mirna_annotation.Rdata")
save(precursor.mirannot, file = "../2_lma96_precursor_mirna_annotation.Rdata")