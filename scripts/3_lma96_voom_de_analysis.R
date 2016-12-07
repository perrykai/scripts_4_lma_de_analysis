#' **Script:** `2_lma96_voom_de_analysis.R`
#' 
#' **Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/4_lma_de_analysis/scripts`
#' 
#' **Date:**  11/08/16; 11/14/16
#' 
#' **Input File Directory:**  
#' 
#' 1. `/mnt/research/pigeqtl/analyses/microRNA/4_lma_de_analysis`
#' 
#' **Input File(s):**
#' 
#' 1. `1_lma96_rounded_mean_mature_mirna_expression.Rdata`
#' 2. `3_lma96_mature_mirna_annotation.Rdata`
#' 
#' **Output File Directory:** 
#' 
#' 1. `/mnt/research/pigeqtl/analyses/microRNA/4_lma_de_analysis/`
#' 
#' **Output File(s):** 
#' 
#' 1. `4_lma96_dge_object_design_mx.Rdata`
#' 2. `5_lma96_voom_de_analysis_results.Rdata`
#' 
#' **Table of contents:**
#'
#' 1. [Objectives](#objectives)
#' 2. [Install libraries](#install-libraries)
#' 3. [Load data](#load-data)
#' 4. [Analysis](#analysis)
#' 5. [Save data](#save-data)

#' ## Objectives
#' The objective of this script is to conduct a differential expression analysis of the 96 animals selected for loin muscle area (LMA) extremes.
#' This analysis will utilize voom and the limma pipeline.
#' 
#' ## Install libraries
library(methods)
library(limma)
library(edgeR)
library(qvalue)

sessionInfo()

#' ## Load data
rm(list=ls())

setwd("/mnt/research/pigeqtl/analyses/microRNA/4_lma_de_analysis/scripts/")
#' Load the miRNA count data, and covariate data
load("../1_lma96_covar_rounded_mean_mature_mirna_expression.Rdata")
#' Load the mature miRNA annotation file
load("../3_lma96_mature_mirna_annotation.Rdata")

ls()

#' ## Analysis
dim(lma96.dfmeanrcround)
head(lma96.dfmeanrcround)

#' ### 1. Create the dge object
#' 
#' First, ensure that the rownames are the same between the counts and annotation files:
if (sum(rownames(lma96.dfmeanrcround)!=rownames(lma96.total.mature.annot2))!=0) stop ("rownames not equal between datasets")

lma96dge<-DGEList(counts=lma96.dfmeanrcround, genes=lma96.total.mature.annot2)
names(lma96dge)
dim(lma96dge$counts)
head(lma96dge$counts)
lma96dge$samples

#' ### 2. Filter the joint expression dataset (dge object)

#' First, eliminate miRNAs whose total expression is 0 among the datasets
sum(rowSums(lma96dge$counts) == 0)

#' So, 79 miRNAs have 0 expression, meaning 332 miRNAs remain

lma96dge<-lma96dge[rowSums(lma96dge$counts)>0,]
dim(lma96dge)

#' Calculate the read counts per million in order to filter miRNAs by normalized expression:
cpm.lma96dge<-cpm(lma96dge)
dim(cpm.lma96dge)
cpm.lma96dge[1:5,1:5]

#' Filter miRNAs with at least 1 cpm in at least 1/4 of the samples (96/4=24)
filtercpm<-rowSums(cpm.lma96dge>=1)>=24
sum(filtercpm)

nrow(cpm.lma96dge) - sum(filtercpm)

#' We are removing 39 miRNA profiles from the analysis
#' 
#' So, keep the miRNA profiles in dge based on those retained in the cpm-filtering step:
#' 
#' This retains the rounded, mean read counts, not the cpm (this will be done later):

lma96dge<-lma96dge[filtercpm,]
names(lma96dge)
lma96dge[1:5,1:5]
dim(lma96dge$counts)

if (sum(colnames(lma96dge)!=colnames(cpm.lma96dge))!=0) stop ("colnames not the same between dge and cpm.lma96dge")

#' Apply the TMM normalization (normalizes for RNA composition by finding a set of scaling factors for the library sizes that minimize the log-fold changes between the samples for the most genes):
#' 
#' The result is the effective library size, which is equal to the product of the original library size and the scaling factor. Effective library size is what is used in downstream analyses.
lma96dge<-calcNormFactors(lma96dge)
lma96dge$samples
hist(lma96dge$samples$norm.factors)
summary(lma96dge$samples$norm.factors)

#' This function (estimateCommonDisp) applies normalization factors, caluclates normalized expression based on robust count of normalized reads.
#' 
#' estimateCommonDisp estimates the quantile-adjusted conditional maximum likelihood (qCML) common dispersion.
#' 
#' This creates a matrix of pseudo-counts, used internally to speed up computation of conditional likelihood used for dispersion estimation & exact tests in the limma pipeline.
#' Pseudo-counts represent the equivalent counts that would have been observed had the library sizes all been equal, assuming the fitted model.
#' DO NOT INTERPRET PSEUDO-COUNTS AS GENERAL-PURPOSE NORMALIZED COUNTS.
lma96dge<-estimateCommonDisp(lma96dge,verbose=TRUE)
lma96dge$common.dispersion

#' Extract the factors needed for the design matrix from the lmapigs covar dataset, these are also used for the MDS plot
sex<-lmapigs$sex
sex<-relevel(sex, ref="M")
head(sex)
status<-lmapigs$Status
status<-relevel(status, ref="L")
head(status)

#' Create MDS plot of the samples, differentiating sex with color (male = Blue, female = Red)
plotMDS(lma96dge, labels=sex,col=ifelse(sex=="M","blue","red"))
legend("bottomright", legend=c("M","F"), pch=20, col=c("blue","red"), ncol=2)

#' Create MDS plot of the samples, differentiating Status with color (High = Green, Low = Red)
plotMDS(lma96dge, labels=status,col=ifelse(status=="L", "red","darkgreen"))
legend("bottomright", legend=c("L","H"), pch=20, col=c("red","darkgreen"), ncol=2)

newcpm.lma96dge<-cpm(lma96dge)
newcpm.lma96dge[1:5,1:5]

#' ### 3. Build the design matrix and apply the voom transformation
#' 
#' First, create factors of sex and status to create the design matrix:
names(lmapigs)
str(lmapigs)

#' Note that "sex" and "Status" are both factors
rownames(lmapigs)<-lmapigs$id
lmapigs[1:5,1:10]

if (sum(colnames(lma96dge)!=rownames(lmapigs))!=0) stop ("rownames not the same between counts and covars datasets")

#' View this as a data.frame to check for correctness:
df.design<-data.frame(counts.id=colnames(lma96dge$counts), lmapigs$id, sex, status)
df.design
#' Just to see, make the sex and status columns numeric and add the combinations of F-H, F-L, M-H, M-L
df.design$sexnum<-as.numeric(df.design$sex)
df.design$statusnum<-as.numeric(df.design$status)
head(df.design)

for (i in 1:nrow(df.design)){
	df.design$sexstatus[i]<-sum(df.design$sexnum[i], df.design$statusnum[i])
}

head(df.design)
table(df.design$sexstatus)


#' Use the model.matrix() command to create the design matrix for this analysis:
design <- model.matrix(~sex+status)
rownames(design)<-colnames(lma96dge$counts)
design

head(lmapigs$Status)
lmapigs$sex[1:10]

#' Apply the voom transformation to the dataset
v <- voom(lma96dge, design=design, plot=TRUE)
names(v)
dim(v$E)
v$E[1:5,1:5]


#' ### 4. Perform the DE analysis using lmFit from limma R package
fit<-lmFit(v, design=design)
names(fit)
head(fit$genes)
#' Residual standard deviations for each gene:
head(fit$sigma)
summary(fit$sigma)

#' Compute moderated t-statistics, moderated F-statistics, and log-odds of DE by empirical Bayes moderation of standard errors towards a common value:
fitmodt<-eBayes(fit)
names(fitmodt)

summary(decideTests(fitmodt))


#' Summarize the differentially expressed genes between the selection statuses using tobTable.
#' The argument "adjust.method = "BH" means using the FDR method to correct for multiple testing.
toptab<-(topTable(fitmodt, coef="statusH", n=Inf, sort="p",adjust.method="BH"))
dim(toptab)
head(toptab)
rownames(toptab)[1:5]

sum(toptab$adj.P.Val <= 0.05)
sum(toptab$adj.P.Val <= 0.01)

hist(toptab$adj.P.Val, ylim=c(0,200))

#' ## Save data
save(lma96dge, design, file="../4_lma96_dge_object_design_mx.Rdata")
save(v, fit, toptab, file="../5_lma96_voom_de_analysis_results.Rdata")