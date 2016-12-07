**Script:** `9_lma96_obtain_target_gene_names.R`

**Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/4_lma_de_analysis/scripts`

**Date:**  12/6/16

**Input File Directory:**  

1. `/mnt/research/pigeqtl/analyses/microRNA/4_lma_de_analysis/`
2. `/mnt/research/pigeqtl/analyses/microRNA/reference_sequences/targetscan_databases`

**Input File(s):** ``

1. `9_targetscan_analysis_results.txt`
2. `Gene_info.txt`

**Output File Directory:** 

`/mnt/research/pigeqtl/analyses/microRNA/4_lma_de_analysis/`

**Output File(s):** 

1. `10_lma96_target_gene_transcript_ids_gene_symbols.txt`
2. `11_lma96_target_transcript_ids.txt`

**Table of contents:**

1. [Objectives](#objectives)
2. [Install libraries](#install-libraries)
3. [Load data](#load-data)
4. [Analysis](#analysis)
5. [Visualize](#visualize)
6. [Save data](#save-data)

## Objectives

The objective of this script is to extract the transcript IDs from the unique targetscan results, and match their appropriate gene symbols for use in the DAVID/KEGG GO analysis

## Install libraries


```r
setwd("/mnt/research/pigeqtl/analyses/microRNA/4_lma_de_analysis/scripts/")
```

## Load data

First, the target gene information


```r
targetgenes<-read.table("../9_targetscan_analysis_results.txt", header=T, sep="\t", colClasses=c(rep("character", 2), rep("numeric", 6), rep("character", 3), rep("numeric", 2)))
dim(targetgenes)
```

```
## [1] 1466326      14
```

```r
head(targetgenes)
```

```
##   a_Gene_ID miRNA_family_ID species_ID MSA_start MSA_end UTR_start UTR_end
## 1    CDR1as     miR-1277-3p       9606      1306    1312       517     523
## 2    CDR1as     miR-1277-3p       9606      2089    2095       736     742
## 3    CDR1as     miR-1277-3p       9606      3170    3194      1242    1247
## 4    CDR1as      miR-129-5p       9606      2104    2167       747     753
## 5    CDR1as      miR-129-5p       9606      2338    2371       856     862
## 6    CDR1as     miR-1306-3p       9606      3654    3659      1424    1429
##   Group_num Site_type miRNA.in.this.species Group_type
## 1         1   7mer-m8                     x    7mer-m8
## 2         2   7mer-m8                     x    7mer-m8
## 3         3      6mer                     x       6mer
## 4         4   7mer-1a                     x    7mer-1a
## 5         5   7mer-1a                     x    7mer-1a
## 6         6      6mer                     x       6mer
##   Species_in_this_group Species_in_this_group_with_this_site_type
## 1                  9606                                        NA
## 2                  9606                                        NA
## 3                  9606                                        NA
## 4                  9606                                        NA
## 5                  9606                                        NA
## 6                  9606                                        NA
##   ORF_overlap
## 1           0
## 2           0
## 3           0
## 4           0
## 5           0
## 6           0
```

Next, the human gene info dataset (contains Gene ID, Gene Symbol, Description, etc)


```r
geneinfo<-read.table("../../reference_sequences/targetscan_databases/Gene_info.txt", header=T, encoding="latin1",quote="",sep="\t", colClasses=c(rep("character",4),rep("numeric", 3)))
dim(geneinfo)
```

```
## [1] 28352     7
```

```r
head(geneinfo)
```

```
##       Transcript.ID            Gene.ID Gene.symbol
## 1            CDR1as             CDR1as      CDR1as
## 2 ENST00000000233.5  ENSG00000004059.6        ARF5
## 3 ENST00000000412.3  ENSG00000003056.3        M6PR
## 4 ENST00000001008.4  ENSG00000004478.5       FKBP4
## 5 ENST00000001146.2  ENSG00000003137.4     CYP26B1
## 6 ENST00000002125.4 ENSG00000003509.11     NDUFAF7
##                                                                              Gene.description
## 1                                                                         circular RNA CDR1as
## 2                                      ADP-ribosylation factor 5 [Source:HGNC Symbol;Acc:658]
## 3               mannose-6-phosphate receptor (cation dependent) [Source:HGNC Symbol;Acc:6752]
## 4                                FK506 binding protein 4, 59kDa [Source:HGNC Symbol;Acc:3720]
## 5       cytochrome P450, family 26, subfamily B, polypeptide 1 [Source:HGNC Symbol;Acc:20581]
## 6 NADH dehydrogenase (ubiquinone) complex I, assembly factor 7 [Source:HGNC Symbol;Acc:28816]
##   Species.ID X3P.seq.tags Representative.transcript.
## 1       9660            1                          1
## 2       9606         2941                          1
## 3       9606           46                          1
## 4       9606          727                          1
## 5       9606         1315                          1
## 6       9606          161                          1
```

## Analysis

First, I want to extract the transcript IDs from the gene targets (column 1)


```r
transcriptids<-unique(targetgenes[,1])
```

How many unique trancript IDs are in the targets?


```r
length(transcriptids)
```

```
## [1] 27430
```

```r
head(transcriptids)
```

```
## [1] "CDR1as"            "ENST00000000233.5" "ENST00000000412.3"
## [4] "ENST00000001008.4" "ENST00000001146.2" "ENST00000002125.4"
```

Use this vector as an index to extract the correct gene name from the gene_info dataset:


```r
genesymbols<-geneinfo[match(transcriptids, geneinfo$Transcript.ID),c("Transcript.ID","Gene.symbol")]
dim(genesymbols)
```

```
## [1] 27430     2
```

```r
head(genesymbols)
```

```
##       Transcript.ID Gene.symbol
## 1            CDR1as      CDR1as
## 2 ENST00000000233.5        ARF5
## 3 ENST00000000412.3        M6PR
## 4 ENST00000001008.4       FKBP4
## 5 ENST00000001146.2     CYP26B1
## 6 ENST00000002125.4     NDUFAF7
```

```r
if(sum(genesymbols$Transcript.ID != transcriptids) != 0) stop ("Transcript ID match index did not work correctly")
```

## Save data


```r
write.table(genesymbols, file="../10_lma96_target_gene_transcript_ids_gene_symbols.txt", quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t")
write.table(genesymbols$Transcript.ID, file="../11_lma96_target_transcript_ids.txt", quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t")
```

