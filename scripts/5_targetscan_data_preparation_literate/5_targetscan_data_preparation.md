**Script:** `5_targetscan_data_preparation.R`

**Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/4_lma_de_analysis/scripts`

**Date:**  11/29/16

**Input File Directory:**  

1. `/mnt/research/pigeqtl/analyses/microRNA/reference_sequences/targetscan_databases`
2. `/mnt/research/pigeqtl/analyses/microRNA/4_lma_de_analysis`

**Input File(s):** 

1. `UTR_sequences_hsa.txt`
2. `miR_Family_Info_hsa.txt`
3. `5_lma96_voom_de_analysis_results.R`

**Output File Directory:** `/mnt/research/pigeqtl/analyses/microRNA/4_lma_de_analysis`

**Output File(s):** 

1. `6_sig_mirna_names.txt`
2. `7_lma96_sig_mirna.txt (output from script 6_lma96_sig_mir_fasta_prep.sh, executed in this script)`
3. `8_lma96_final_sig_hsa_mirna.txt`

**Table of contents:**

1. [Objectives](#objectives)
2. [Install libraries](#install-libraries)
3. [Load data](#load-data)
4. [Analysis](#analysis)
5. [Visualize](#visualize)
6. [Save data](#save-data)

## Objectives

The objective of this script is to prepare the miRNA datasets for use in the targetscan target prediction analysis

The list of significant miRNAs identified in the DE analysis will be extracted and used as an index to subset the hsa reference datasets.

This will allow the human-equivalent miRNAs to be queried against the human gene list for target prediction. 

Another point to consider would be to limit the UTR data to only those genes expressed in the 174 F2 pigs, using a gene list provided by Deborah.

## Install libraries
## Load data


```r
setwd("/mnt/research/pigeqtl/analyses/microRNA/4_lma_de_analysis/scripts")
```

Load the results from the DE Analysis


```r
load("../5_lma96_voom_de_analysis_results.Rdata")
```

Load the human miRNAs


```r
hsamirna<-read.table("../../reference_sequences/targetscan_databases/miR_Family_Info_hsa.txt")
colnames(hsamirna)<-c("mirna_name", "seed_seq", "species_id")
head(hsamirna)
```

```
##   mirna_name seed_seq species_id
## 1   miR-9500  AGGGAAG       9606
## 2   miR-8485  ACACACA       9606
## 3   miR-8088  CUCGGUA       9606
## 4   miR-8087  AAGACUU       9606
## 5   miR-8086  GCUAGUC       9606
## 6   miR-8084  AAUACUA       9606
```

```r
hsamirna$mirna_name<-as.character(hsamirna$mirna_name)
```

## Analysis

Isolate the names of the miRNAs significant in the DE analysis


```r
sum(toptab$adj.P.Val < 0.05)
```

```
## [1] 148
```

```r
sigmir<-toptab[toptab$adj.P.Val < 0.05,]
dim(sigmir)
```

```
## [1] 148  15
```

```r
pigmirnm<-as.character(sigmir$Name)
length(pigmirnm)
```

```
## [1] 148
```

```r
pigmir<-sub("ssc-", "", pigmirnm)
head(pigmir)
```

```
## [1] "miR-99b"     "miR-191"     "miR-136"     "miR-376a-3p" "miR-19b"    
## [6] "miR-19a"
```

The following allows the significant pig miRNAs to be extracted correctly from the full dataset to compare seed sequences to human miR


```r
pigmirnm2<-paste(pigmirnm, '$', sep="")
```

Write this to a .txt file right now, then execute the next script from this script using that file


```r
write.table(pigmirnm2, file="../6_sig_mirna_names.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)
```

Change the permissions on the next script:


```r
system("chmod 774 6_lma96_sig_mir_fasta_prep.sh")
```

Execute the next script so the analysis can continue


```r
system("./6_lma96_sig_mir_fasta_prep.sh", wait=FALSE)
```

Need to remove the "3p", "5p" information to match between this dataset and the human dataset
(this is ugly, but it works - ask if there's a better way to do this, haha)


```r
pigsub<-strsplit(pigmir, "-")
pigsub<-do.call(rbind, pigsub)[,1:2]
```

```
## Warning in (function (..., deparse.level = 1) : number of columns of result
## is not a multiple of vector length (arg 4)
```

```r
head(pigsub)
```

```
##      [,1]  [,2]  
## [1,] "miR" "99b" 
## [2,] "miR" "191" 
## [3,] "miR" "136" 
## [4,] "miR" "376a"
## [5,] "miR" "19b" 
## [6,] "miR" "19a"
```

```r
pigsub<-paste(pigsub[,1], pigsub[,2], sep="-")
head(pigsub)
```

```
## [1] "miR-99b"  "miR-191"  "miR-136"  "miR-376a" "miR-19b"  "miR-19a"
```

```r
length(pigsub)
```

```
## [1] 148
```

How many unique miRNAs are in the dataset?


```r
length(unique(pigsub))
```

```
## [1] 131
```

Use these names to find matches against the human miRNA database


```r
y<-list()
x<-list()
hsasub<-data.frame()
for(i in pigsub){
	y[[i]]<-grep(i, hsamirna$mirna_name)
	x[[i]]<-hsamirna[y[[i]],]
}
length(y)
```

```
## [1] 131
```

```r
length(x)
```

```
## [1] 131
```

```r
xdf<-do.call(rbind, x)
head(xdf)
```

```
##                       mirna_name seed_seq species_id
## miR-191.1354         miR-1915-5p  CCUUGCC       9606
## miR-191.1355 miR-1915-3p/6764-5p  CCCAGGG       9606
## miR-191.1356 miR-1915-3p/6764-5p  CCCAGGG       9606
## miR-191.1357         miR-1911-5p  GAGUACC       9606
## miR-191.1358         miR-1914-5p  CCUGUGC       9606
## miR-191.1359         miR-1910-5p  CAGUCCU       9606
```

```r
dim(xdf)
```

```
## [1] 880   3
```

How many miRNAs had a match in this step? (i.e., how many am I gonna have to go and find myself?)


```r
length(unique(gsub("\\..*","",rownames(xdf))))
```

```
## [1] 103
```

Instead of using the names to match the human miRNAs, use the seed sequences!

Load the known pig miRNA sequences corresponding to the significantly DE miRNAs


```r
pigseq<-read.table("../7_lma96_sig_mirna.txt", col.names=c("name","seq"))
head(pigseq)
```

```
##               name                     seq
## 1    >ssc-miR-125b  UCCCUGAGACCCUAACUUGUGA
## 2 >ssc-miR-148a-5p  AAAGUUCUGAGACACUCCGACU
## 3     >ssc-miR-15b  UAGCAGCACAUCAUGGUUUACA
## 4     >ssc-miR-184  UGGACGGAGAACUGAUAAGGGU
## 5     >ssc-miR-19a UGUGCAAAUCUAUGCAAAACUGA
## 6     >ssc-miR-20a  UAAAGUGCUUAUAGUGCAGGUA
```

```r
pigseq$name<-sub(">", "", pigseq$name)
```

Extract nt 2-8 (the miRNA seed sequence)


```r
pigseq$seedseq<-substr(pigseq$seq, 2, 8)
head(pigseq)
```

```
##              name                     seq seedseq
## 1    ssc-miR-125b  UCCCUGAGACCCUAACUUGUGA CCCUGAG
## 2 ssc-miR-148a-5p  AAAGUUCUGAGACACUCCGACU AAGUUCU
## 3     ssc-miR-15b  UAGCAGCACAUCAUGGUUUACA AGCAGCA
## 4     ssc-miR-184  UGGACGGAGAACUGAUAAGGGU GGACGGA
## 5     ssc-miR-19a UGUGCAAAUCUAUGCAAAACUGA GUGCAAA
## 6     ssc-miR-20a  UAAAGUGCUUAUAGUGCAGGUA AAAGUGC
```

Match the seed sequences between the pig seq and the human seq, returning the miRNA name of the matching human seq


```r
pigseq$matchedseed<-xdf[match(pigseq$seedseq, xdf$seed_seq),"mirna_name"]
head(pigseq)
```

```
##              name                     seq seedseq
## 1    ssc-miR-125b  UCCCUGAGACCCUAACUUGUGA CCCUGAG
## 2 ssc-miR-148a-5p  AAAGUUCUGAGACACUCCGACU AAGUUCU
## 3     ssc-miR-15b  UAGCAGCACAUCAUGGUUUACA AGCAGCA
## 4     ssc-miR-184  UGGACGGAGAACUGAUAAGGGU GGACGGA
## 5     ssc-miR-19a UGUGCAAAUCUAUGCAAAACUGA GUGCAAA
## 6     ssc-miR-20a  UAAAGUGCUUAUAGUGCAGGUA AAAGUGC
##                            matchedseed
## 1                           miR-125-5p
## 2                          miR-148a-5p
## 3 miR-15-5p/16-5p/195-5p/424-5p/497-5p
## 4                              miR-184
## 5                            miR-19-3p
## 6  miR-17-5p/20-5p/93-5p/106-5p/519-3p
```

```r
head(pigseq[,c(1,4)])
```

```
##              name                          matchedseed
## 1    ssc-miR-125b                           miR-125-5p
## 2 ssc-miR-148a-5p                          miR-148a-5p
## 3     ssc-miR-15b miR-15-5p/16-5p/195-5p/424-5p/497-5p
## 4     ssc-miR-184                              miR-184
## 5     ssc-miR-19a                            miR-19-3p
## 6     ssc-miR-20a  miR-17-5p/20-5p/93-5p/106-5p/519-3p
```

```r
sum(is.na(pigseq$matchedseed))
```

```
## [1] 47
```

This way matches the seed sequences between the full dataset of human miRNAs, not the filtered


```r
pigseq$matchedseed2<-hsamirna[match(pigseq$seedseq, hsamirna$seed_seq), "mirna_name"]
head(pigseq[,c(1,5)])
```

```
##              name                         matchedseed2
## 1    ssc-miR-125b                           miR-125-5p
## 2 ssc-miR-148a-5p                          miR-148a-5p
## 3     ssc-miR-15b miR-15-5p/16-5p/195-5p/424-5p/497-5p
## 4     ssc-miR-184                              miR-184
## 5     ssc-miR-19a                            miR-19-3p
## 6     ssc-miR-20a  miR-17-5p/20-5p/93-5p/106-5p/519-3p
```

```r
sum(is.na(pigseq$matchedseed2))
```

```
## [1] 20
```

```r
sum(pigseq$matchedseed2 %in% pigseq$matchedseed)
```

```
## [1] 121
```

```r
dim(pigseq)
```

```
## [1] 148   5
```

```r
nrow(pigseq) - sum(pigseq$matchedseed2 %in% pigseq$matchedseed)
```

```
## [1] 27
```

```r
sum(is.na(pigseq$matchedseed)) - sum(is.na(pigseq$matchedseed2))
```

```
## [1] 27
```

```r
head(pigseq)
```

```
##              name                     seq seedseq
## 1    ssc-miR-125b  UCCCUGAGACCCUAACUUGUGA CCCUGAG
## 2 ssc-miR-148a-5p  AAAGUUCUGAGACACUCCGACU AAGUUCU
## 3     ssc-miR-15b  UAGCAGCACAUCAUGGUUUACA AGCAGCA
## 4     ssc-miR-184  UGGACGGAGAACUGAUAAGGGU GGACGGA
## 5     ssc-miR-19a UGUGCAAAUCUAUGCAAAACUGA GUGCAAA
## 6     ssc-miR-20a  UAAAGUGCUUAUAGUGCAGGUA AAAGUGC
##                            matchedseed
## 1                           miR-125-5p
## 2                          miR-148a-5p
## 3 miR-15-5p/16-5p/195-5p/424-5p/497-5p
## 4                              miR-184
## 5                            miR-19-3p
## 6  miR-17-5p/20-5p/93-5p/106-5p/519-3p
##                           matchedseed2
## 1                           miR-125-5p
## 2                          miR-148a-5p
## 3 miR-15-5p/16-5p/195-5p/424-5p/497-5p
## 4                              miR-184
## 5                            miR-19-3p
## 6  miR-17-5p/20-5p/93-5p/106-5p/519-3p
```

The next step is to use the name of the human miRNA sequences to subset the larger database for use with the targetscan algorithm


```r
head(hsamirna)
```

```
##   mirna_name seed_seq species_id
## 1   miR-9500  AGGGAAG       9606
## 2   miR-8485  ACACACA       9606
## 3   miR-8088  CUCGGUA       9606
## 4   miR-8087  AAGACUU       9606
## 5   miR-8086  GCUAGUC       9606
## 6   miR-8084  AAUACUA       9606
```

```r
finalhsamir<-hsamirna[match(pigseq$matchedseed2, hsamirna$mirna_name),]
finalhsamir<-finalhsamir[complete.cases(finalhsamir),]
dim(finalhsamir)
```

```
## [1] 128   3
```

```r
head(finalhsamir)
```

```
##                                mirna_name seed_seq species_id
## 2417                           miR-125-5p  CCCUGAG       9606
## 2341                          miR-148a-5p  AAGUUCU       9606
## 2560 miR-15-5p/16-5p/195-5p/424-5p/497-5p  AGCAGCA       9606
## 2303                              miR-184  GGACGGA       9606
## 2539                            miR-19-3p  GUGCAAA       9606
## 2546  miR-17-5p/20-5p/93-5p/106-5p/519-3p  AAAGUGC       9606
```

## Save data


```r
write.table(finalhsamir, file="../8_lma96_final_sig_hsa_mirna.txt", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
```

