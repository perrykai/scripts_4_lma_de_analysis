**Script:** `1_extract_covar_expression_data.R`

**Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/4_lma_de_analysis/scripts`

**Date:**  11/8/16

**Input File Directory:**  

1. `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/3_build_dge_object_for_eqtl/`
2. `/mnt/research/pigeqtl/analyses/microRNA/1_preprocess_fastq_files/10_mirdeep2_core_quantify_predict_output/`
3. `/mnt/research/pigeqtl/analyses/microRNA/1_preprocess_fastq_files/8_collapsed_fasta_output_expression_matrix/`

**Input File(s):** 

1. `4_msuprp_mirna_gpdata_pheno_counts.Rdata`
2. `miRNAs_expressed_all_samples_21_01_2016_t_20_01_10.csv`
3. `config_for_mapper.txt`

**Output File Directory:** `/mnt/research/pigeqtl/analyses/microRNA/4_lma_de_analysis/`

**Output File(s):** `1_lma96_rounded_mean_mature_mirna_expression.Rdata`

**Table of contents:**

1. [Objectives](#objectives)
2. [Install libraries](#install-libraries)
3. [Load data](#load-data)
4. [Analysis](#analysis)
5. [Visualize](#visualize)
6. [Save data](#save-data)

## Objectives

The objective of this script is to extract the miRNA expression data pertaining to the 96 animals that will be used in the miRNA DE analysis for LMA extremes.

## Install libraries


```r
rm(list=ls())
library(plyr)
```

## Load data


```r
setwd("/mnt/research/pigeqtl/analyses/microRNA/4_lma_de_analysis/scripts")

load("../../2_mirna_characterization_expression/3_build_dge_object_for_eqtl/4_msuprp_mirna_gpdata_pheno_counts.Rdata")
```

Read in the csv of the miRNA expression profiles, using check.names default to correct column names with non-ASCII characters


```r
rc<-read.csv("../../1_preprocess_fastq_files/10_mirdeep2_core_quantify_predict_output/miRNAs_expressed_all_samples_21_01_2016_t_20_01_10.csv", sep = "\t", header = TRUE, row.names=NULL)
```

Read in the config file, maintaining the characters in the 3-digit code names


```r
configfile<-read.table("../../1_preprocess_fastq_files/8_collapsed_fasta_output_expression_matrix/config_for_mapper.txt", header = FALSE, sep = " ", row.names=NULL, colClasses = c('character','character'))
head(configfile)
```

```
##                V1  V2
## 1 1034_express.fa 001
## 2 1036_express.fa 002
## 3 1041_express.fa 003
## 4 1049_express.fa 004
## 5 1058_express.fa 005
## 6 1060_express.fa 006
```

```r
dim(configfile)
```

```
## [1] 174   2
```

Remove the "_express.fa" from each file name to leave the pig id:


```r
configfile$V1<-gsub("_express.fa", "", configfile$V1)
```

Make filenames more informative:


```r
colnames(configfile)<-c("pigid","code")
head(configfile)
```

```
##   pigid code
## 1  1034  001
## 2  1036  002
## 3  1041  003
## 4  1049  004
## 5  1058  005
## 6  1060  006
```

```r
colnames(configfile)
```

```
## [1] "pigid" "code"
```

```r
ls()
```

```
## [1] "configfile"                 "final_MSUPRP_miRNA"        
## [3] "rc"                         "summary_final_MSUPRP_miRNA"
```

## Analysis

First, extract the IDs from the animals meeting the selection criteria of "lma", loin muscle area.


```r
names(final_MSUPRP_miRNA)
```

```
## [1] "covar"       "pheno"       "geno"        "map"         "pedigree"   
## [6] "phenoCovars" "info"
```

```r
head(final_MSUPRP_miRNA$covar)
```

```
##      id phenotyped genotyped family sex litter wt_10wk bf10_22wk lma_22wk
## 26 1034       TRUE      TRUE     NA   F      5   30.84     16.51    47.55
## 28 1036       TRUE      TRUE     NA   F      5   29.48     18.80    31.16
## 33 1041       TRUE      TRUE     NA   M      5   32.66     21.08    41.94
## 41 1049       TRUE      TRUE     NA   M      5   29.03     21.59    31.94
## 50 1058       TRUE      TRUE     NA   F     10   24.95     10.16    39.10
## 52 1060       TRUE      TRUE     NA   F     10   20.87     17.02    31.03
##    slgdt_cd age_slg car_wt Status selcrit microarray_Dye microarray_file
## 26        2     160  85.49      H     lma            cy3     slide33.gpr
## 28        2     160  71.66      L     lma            cy5     slide33.gpr
## 33        2     160  80.27      H     lma            cy5     slide34.gpr
## 41        2     160  81.86      L     lma            cy3     slide34.gpr
## 50        2     159  72.11      H     lma            cy5     slide52.gpr
## 52        2     159  70.75      L     lma            cy3     slide52.gpr
##    Color Hairden earset Spots Underbelly face_color line_origin perc_duroc
## 26  <NA>    <NA>   <NA>  <NA>       <NA>         NA          NA  0.5051327
## 28  <NA>    <NA>   <NA>  <NA>       <NA>         NA          NA  0.5471133
## 33  <NA>    <NA>   <NA>  <NA>       <NA>         NA          NA  0.4721328
## 41  <NA>    <NA>   <NA>  <NA>       <NA>         NA          NA  0.4565887
## 50  <NA>    <NA>   <NA>  <NA>       <NA>         NA          NA  0.4336773
## 52  <NA>    <NA>   <NA>  <NA>       <NA>         NA          NA  0.5354381
##    growth_group
## 26        lma-H
## 28        lma-L
## 33        lma-H
## 41        lma-L
## 50        lma-H
## 52        lma-L
```

```r
table(final_MSUPRP_miRNA$covar$selcrit)
```

```
## 
##  bf lma 
##  78  96
```

```r
final_MSUPRP_miRNA$covar$selcrit
```

```
##   [1] lma lma lma lma lma lma lma lma lma lma lma lma lma lma lma lma lma
##  [18] lma lma lma lma lma lma lma lma lma lma lma lma lma lma lma lma lma
##  [35] lma lma lma lma lma lma lma lma lma lma lma lma lma lma lma lma lma
##  [52] lma lma lma lma lma lma lma lma lma lma lma lma lma lma lma lma lma
##  [69] lma lma lma lma lma lma lma lma lma lma lma lma lma lma lma lma lma
##  [86] lma lma lma lma lma lma lma lma lma lma lma bf  bf  bf  bf  bf  bf 
## [103] bf  bf  bf  bf  bf  bf  bf  bf  bf  bf  bf  bf  bf  bf  bf  bf  bf 
## [120] bf  bf  bf  bf  bf  bf  bf  bf  bf  bf  bf  bf  bf  bf  bf  bf  bf 
## [137] bf  bf  bf  bf  bf  bf  bf  bf  bf  bf  bf  bf  bf  bf  bf  bf  bf 
## [154] bf  bf  bf  bf  bf  bf  bf  bf  bf  bf  bf  bf  bf  bf  bf  bf  bf 
## [171] bf  bf  bf  bf 
## Levels: bf lma
```

```r
lmapigs<-final_MSUPRP_miRNA$covar[final_MSUPRP_miRNA$covar$selcrit == "lma",]
dim(lmapigs)
```

```
## [1] 96 25
```

```r
table(lmapigs$Status)
```

```
## 
##  H  L 
## 48 48
```

```r
table(lmapigs$sex)
```

```
## 
##  F  M 
## 48 48
```

```r
ids<-lmapigs$id

head(rc[,1:5])
```

```
##         X.miRNA read_count    precursor   total  X001
## 1    ssc-let-7a    5066155 ssc-let-7a-1 5066155 48170
## 2    ssc-let-7a    5055934 ssc-let-7a-2 5055934 48095
## 3    ssc-let-7c    3238382   ssc-let-7c 3238382 32745
## 4 ssc-let-7d-5p     490052   ssc-let-7d  490052  4925
## 5 ssc-let-7d-3p      51589   ssc-let-7d   51589   381
## 6    ssc-let-7e     274145   ssc-let-7e  274145  2811
```

Set the name of the first column to "miRNA":


```r
colnames(rc)[[1]]<-"miRNA"
```

Remove the "X" character from the beginning of each column:


```r
colnames(rc)<-gsub("X", "", colnames(rc))
colnames(rc)[1:30]
```

```
##  [1] "miRNA"      "read_count" "precursor"  "total"      "001"       
##  [6] "002"        "003"        "004"        "005"        "006"       
## [11] "007"        "008"        "009"        "010"        "011"       
## [16] "012"        "013"        "014"        "015"        "016"       
## [21] "017"        "018"        "019"        "020"        "021"       
## [26] "022"        "023"        "024"        "025"        "026"
```

### 1. Extract the columns of mature read counts for each miRNA for each animal


```r
lmaquant<-rc[,c(1,5:178)]
colnames(lmaquant)
```

```
##   [1] "miRNA" "001"   "002"   "003"   "004"   "005"   "006"   "007"  
##   [9] "008"   "009"   "010"   "011"   "012"   "013"   "014"   "015"  
##  [17] "016"   "017"   "018"   "019"   "020"   "021"   "022"   "023"  
##  [25] "024"   "025"   "026"   "027"   "028"   "029"   "030"   "031"  
##  [33] "032"   "033"   "034"   "035"   "036"   "037"   "038"   "039"  
##  [41] "040"   "041"   "042"   "043"   "044"   "045"   "046"   "047"  
##  [49] "048"   "049"   "050"   "051"   "052"   "053"   "054"   "055"  
##  [57] "056"   "057"   "058"   "059"   "060"   "061"   "062"   "063"  
##  [65] "064"   "065"   "066"   "067"   "068"   "069"   "070"   "071"  
##  [73] "072"   "073"   "074"   "075"   "076"   "077"   "078"   "079"  
##  [81] "080"   "081"   "082"   "083"   "084"   "085"   "086"   "087"  
##  [89] "088"   "089"   "090"   "091"   "092"   "093"   "094"   "095"  
##  [97] "096"   "097"   "098"   "099"   "100"   "101"   "102"   "103"  
## [105] "104"   "105"   "106"   "107"   "108"   "109"   "110"   "111"  
## [113] "112"   "113"   "114"   "115"   "116"   "117"   "118"   "119"  
## [121] "120"   "121"   "122"   "123"   "124"   "125"   "126"   "127"  
## [129] "128"   "129"   "130"   "131"   "132"   "133"   "134"   "135"  
## [137] "136"   "137"   "138"   "139"   "140"   "141"   "142"   "143"  
## [145] "144"   "145"   "146"   "147"   "148"   "149"   "150"   "151"  
## [153] "152"   "153"   "154"   "155"   "156"   "157"   "158"   "159"  
## [161] "160"   "161"   "162"   "163"   "164"   "165"   "166"   "167"  
## [169] "168"   "169"   "170"   "171"   "172"   "173"   "174"
```

```r
dim(lmaquant)
```

```
## [1] 492 175
```

```r
head(lmaquant[1:8])
```

```
##           miRNA   001   002   003   004   005   006   007
## 1    ssc-let-7a 48170 23457 28447 29899 40817 29774 38804
## 2    ssc-let-7a 48095 23397 28448 29821 40698 29725 38794
## 3    ssc-let-7c 32745 14987 18144 18681 34313 15294 28022
## 4 ssc-let-7d-5p  4925  1938  2511  3076  3472  3276  3705
## 5 ssc-let-7d-3p   381   192   198   269   778   239   774
## 6    ssc-let-7e  2811  1302  1463  1690  2512  1410  2229
```

Take a subset of this data.frame for testing:


```r
test<-lmaquant[1:20,1:8]

test[1:10,]
```

```
##            miRNA   001   002   003   004   005   006   007
## 1     ssc-let-7a 48170 23457 28447 29899 40817 29774 38804
## 2     ssc-let-7a 48095 23397 28448 29821 40698 29725 38794
## 3     ssc-let-7c 32745 14987 18144 18681 34313 15294 28022
## 4  ssc-let-7d-5p  4925  1938  2511  3076  3472  3276  3705
## 5  ssc-let-7d-3p   381   192   198   269   778   239   774
## 6     ssc-let-7e  2811  1302  1463  1690  2512  1410  2229
## 7     ssc-let-7f 35638 16475 20241 23079 16391 24685 19165
## 8     ssc-let-7f 35227 16369 20082 22892 16634 24266 19185
## 9     ssc-let-7g 16700  8010  9774 10902 11420 11251 12589
## 10    ssc-let-7i 10966  4496  5306  6397  6161  5474  5958
```

### 2. Use the 'by' function to apply the function colMeans to the entire data frame of read counts:
(What this will do is go down the columns looking at the index of grouped miRNA names and take the average of the read counts for that miRNA)
The result of this will be a list containing the average read counts for each miRNA for each animal.

Example: by(data, index, function)


```r
bytst<-by(test[,2:ncol(test)], test[,1], colMeans)
bytst[1:25]
```

```
## $`ssc-let-7a`
##     001     002     003     004     005     006     007 
## 48132.5 23427.0 28447.5 29860.0 40757.5 29749.5 38799.0 
## 
## $`ssc-let-7c`
##   001   002   003   004   005   006   007 
## 32745 14987 18144 18681 34313 15294 28022 
## 
## $`ssc-let-7d-3p`
## 001 002 003 004 005 006 007 
## 381 192 198 269 778 239 774 
## 
## $`ssc-let-7d-5p`
##  001  002  003  004  005  006  007 
## 4925 1938 2511 3076 3472 3276 3705 
## 
## $`ssc-let-7e`
##  001  002  003  004  005  006  007 
## 2811 1302 1463 1690 2512 1410 2229 
## 
## $`ssc-let-7f`
##     001     002     003     004     005     006     007 
## 35432.5 16422.0 20161.5 22985.5 16512.5 24475.5 19175.0 
## 
## $`ssc-let-7g`
##   001   002   003   004   005   006   007 
## 16700  8010  9774 10902 11420 11251 12589 
## 
## $`ssc-let-7i`
##   001   002   003   004   005   006   007 
## 10966  4496  5306  6397  6161  5474  5958 
## 
## $`ssc-miR-1`
##    001    002    003    004    005    006    007 
## 307458 228189 261470 283840 109839 282455 117277 
## 
## $`ssc-miR-100`
##   001   002   003   004   005   006   007 
##  9539  3194  3066  3069 17474  2405 12515 
## 
## $`ssc-miR-101`
##     001     002     003     004     005     006     007 
## 12024.0  7469.5  9833.0 10507.5  4269.0  5779.0  5073.5 
## 
## $`ssc-miR-103`
##    001    002    003    004    005    006    007 
## 2244.5 1182.0 1309.0 1822.5 1881.5 1044.0 1742.0 
## 
## $`ssc-miR-105-1`
## 001 002 003 004 005 006 007 
##   3   1   0   2   2   0   0 
## 
## $`ssc-miR-105-2`
## 001 002 003 004 005 006 007 
##   1   0   0   0   0   0   0 
## 
## $`ssc-miR-106a`
## 001 002 003 004 005 006 007 
##  39  17  26  45  27  25  22 
## 
## $`ssc-miR-107`
## 001 002 003 004 005 006 007 
## 460 213 256 313 357 136 282 
## 
## $`ssc-miR-10a-3p`
## NULL
## 
## $`ssc-miR-10a-5p`
## NULL
## 
## $`ssc-miR-10b`
## NULL
## 
## $`ssc-miR-122`
## NULL
## 
## $`ssc-miR-1224`
## NULL
## 
## $`ssc-miR-1249`
## NULL
## 
## $`ssc-miR-124a`
## NULL
## 
## $`ssc-miR-125a`
## NULL
## 
## $`ssc-miR-125b`
## NULL
```

Notice here that the rest of the miRNA names remain since the miRNA name is a factor, but since there is no data for them they are filled with NULL.

Apply the by function to the full dataframe:


```r
meanrc<-by(lmaquant[,2:ncol(lmaquant)], lmaquant[,1], colMeans)
```

This should be 411 (the number of mature pig miRNAs in miRBase), meaning we have one expression profile for each mature miRNA:


```r
length(meanrc)
```

```
## [1] 411
```

```r
meanrc[1:3]
```

```
## $`ssc-let-7a`
##     001     002     003     004     005     006     007     008     009 
## 48132.5 23427.0 28447.5 29860.0 40757.5 29749.5 38799.0 31459.5 29857.0 
##     010     011     012     013     014     015     016     017     018 
## 21333.0 35977.0 21304.0 43336.5 46963.0 41826.5 41174.0 36678.5 45540.0 
##     019     020     021     022     023     024     025     026     027 
## 59997.0 34675.0 36866.0 34397.5 42057.0 49982.0 35410.5  2366.0 25973.5 
##     028     029     030     031     032     033     034     035     036 
## 41518.5 42680.0 30864.0 44394.0 43819.0 38598.5 50949.5 46225.0 47411.0 
##     037     038     039     040     041     042     043     044     045 
## 35912.5 37383.0 24522.0 35894.5 15876.5 42837.0 36219.0 30672.0 28829.0 
##     046     047     048     049     050     051     052     053     054 
## 28805.0 24042.0 30291.5 25858.0 31424.0 10909.0 33168.0 23342.0 28218.0 
##     055     056     057     058     059     060     061     062     063 
## 17887.5 15043.5 16239.5 22696.0 13045.5 28524.0 17197.0 32303.5 28680.0 
##     064     065     066     067     068     069     070     071     072 
## 37439.0 47521.5 42747.0 41810.5 20253.0 28069.0 47363.0 29313.0 26772.0 
##     073     074     075     076     077     078     079     080     081 
## 35965.5 27357.5 42022.0 44591.0 48267.5 35515.5 30704.5 37564.5 38083.5 
##     082     083     084     085     086     087     088     089     090 
## 36868.5 40084.0 35862.5 23067.5 37179.0 38866.5 29532.5 28213.0 39304.5 
##     091     092     093     094     095     096     097     098     099 
## 37340.5 33901.5 25400.0 28632.0 26621.0 23326.0 27332.5 29679.5 24334.5 
##     100     101     102     103     104     105     106     107     108 
## 25708.0 29901.0 30335.0 23033.5 23954.0 17311.0 19595.5 28878.5 22962.5 
##     109     110     111     112     113     114     115     116     117 
## 18231.0 25629.5 17310.5 14579.0 13451.5 14827.5 35094.0 22924.0 22506.5 
##     118     119     120     121     122     123     124     125     126 
## 21439.0 19489.5 22782.0 24263.5 21087.0 23661.5 25093.0 25866.0 20312.5 
##     127     128     129     130     131     132     133     134     135 
## 23853.0 28593.0 28171.0 27868.5 23409.5 21416.0 18956.5 29712.0 20679.0 
##     136     137     138     139     140     141     142     143     144 
## 22174.0 23523.5 26256.5 23873.5 23470.0 16744.5 23586.5 32558.0 23146.5 
##     145     146     147     148     149     150     151     152     153 
## 13914.0 15831.0 30359.5 21163.5 17806.5 32457.5 28573.0 26978.5 15400.5 
##     154     155     156     157     158     159     160     161     162 
## 20649.0 20217.0 34496.0 27724.0 32625.0 27000.0 31787.0 25896.0 24638.0 
##     163     164     165     166     167     168     169     170     171 
## 22857.0 23807.0 19823.5 28145.5 34188.5 29473.5 16803.0 24594.5 21326.0 
##     172     173     174 
## 19700.5 20084.0 18977.0 
## 
## $`ssc-let-7c`
##   001   002   003   004   005   006   007   008   009   010   011   012 
## 32745 14987 18144 18681 34313 15294 28022 19512 16741 12522 20772 12449 
##   013   014   015   016   017   018   019   020   021   022   023   024 
## 26354 28968 28788 29017 29022 31957 45232 23212 23731 21229 29962 37824 
##   025   026   027   028   029   030   031   032   033   034   035   036 
## 23684  3077 15389 28463 30503 17962 31221 31310 25291 37902 32444 32430 
##   037   038   039   040   041   042   043   044   045   046   047   048 
## 24015 25399 14477 25063 14181 31924 24089 20474 21753 18148 14466 18737 
##   049   050   051   052   053   054   055   056   057   058   059   060 
## 16355 19779  5973 20260 14482 18147 10798  8526  9964 15666  7201 17704 
##   061   062   063   064   065   066   067   068   069   070   071   072 
##  9342 24491 20807 26015 28812 26738 26126 17415 16928 34465 17004 15809 
##   073   074   075   076   077   078   079   080   081   082   083   084 
## 22883 16277 27810 31228 34068 25721 18956 24904 24512 22762 25993 24232 
##   085   086   087   088   089   090   091   092   093   094   095   096 
## 14980 19671 24475 18429 18530 25909 23004 20516 14806 15282 13519 13682 
##   097   098   099   100   101   102   103   104   105   106   107   108 
## 15386 17262 13800 16006 18661 20238 13522 14096 10464 10629 16443 14089 
##   109   110   111   112   113   114   115   116   117   118   119   120 
## 10157 15627  8287  7647  6847  8201 22353 13792 12973 12208 11024 13482 
##   121   122   123   124   125   126   127   128   129   130   131   132 
## 12707 12936 14476 15825 14335 11533 15121 18384 18005 15034 13522 12939 
##   133   134   135   136   137   138   139   140   141   142   143   144 
## 10933 21979 13147 13539 14349 18145 15222 14097 10471 15412 21981 13851 
##   145   146   147   148   149   150   151   152   153   154   155   156 
##  8184  8413 19455 13359  9886 19750 17841 16159  8110 12113 12460 24578 
##   157   158   159   160   161   162   163   164   165   166   167   168 
## 18455 22931 16695 20916 14315 14228 13955 13703 10785 16482 23031 17416 
##   169   170   171   172   173   174 
##  9464 14097 11644 11737 13957 10290 
## 
## $`ssc-let-7d-3p`
##  001  002  003  004  005  006  007  008  009  010  011  012  013  014  015 
##  381  192  198  269  778  239  774  355  243  166  440  182  389  503  514 
##  016  017  018  019  020  021  022  023  024  025  026  027  028  029  030 
##  449  710  460 1071  411  356  330  319  750  449  525  339  378  393  256 
##  031  032  033  034  035  036  037  038  039  040  041  042  043  044  045 
##  438  407  352  498  569  527  269  326  199  265 1033  365  261  268  652 
##  046  047  048  049  050  051  052  053  054  055  056  057  058  059  060 
##  195  209  203  183  236   84  170  167  212  149   84  101  186   69  147 
##  061  062  063  064  065  066  067  068  069  070  071  072  073  074  075 
##  116  566  308  379  383  536  323  979  264  506  254  161  310  180  487 
##  076  077  078  079  080  081  082  083  084  085  086  087  088  089  090 
##  507  460  492  252  345  369  334  502  292  214  271  306  317  249  417 
##  091  092  093  094  095  096  097  098  099  100  101  102  103  104  105 
##  405  265  295  170  293  245  327  212  246  302  278  293  183  225  130 
##  106  107  108  109  110  111  112  113  114  115  116  117  118  119  120 
##  192  245  211  162  292  150   96   89  113  381  281  224  227  174  152 
##  121  122  123  124  125  126  127  128  129  130  131  132  133  134  135 
##  169  245  276  233  360  255  172  249  224  241  258  211  206  239  184 
##  136  137  138  139  140  141  142  143  144  145  146  147  148  149  150 
##  236  202  233  178  212  180  195  255  219  242   95  254  221  220  439 
##  151  152  153  154  155  156  157  158  159  160  161  162  163  164  165 
##  263  205   99  204  157  333  196  375  265  250  239  203  184  168  178 
##  166  167  168  169  170  171  172  173  174 
##  268  183  196  195  245  204  141  202  153
```


### 3. Transform the list output back into a data.frame using the plyr package to prepare for gblup function:
Example: ldply(.data, .fun, .id)

id = name of the index column (used if data is a named list). Pass NULL to avoid creation
     of the index column. For compatibility, omit this argument or pass "NA" to avoid converting the index column
     to a factor; in this case, ".id" is used as column name.


```r
lma.dfmeanrc<-ldply(meanrc, fun=NULL, id=names(meanrc))

head(lma.dfmeanrc[1:8])
```

```
##             .id     001   002     003     004     005     006   007
## 1    ssc-let-7a 48132.5 23427 28447.5 29860.0 40757.5 29749.5 38799
## 2    ssc-let-7c 32745.0 14987 18144.0 18681.0 34313.0 15294.0 28022
## 3 ssc-let-7d-3p   381.0   192   198.0   269.0   778.0   239.0   774
## 4 ssc-let-7d-5p  4925.0  1938  2511.0  3076.0  3472.0  3276.0  3705
## 5    ssc-let-7e  2811.0  1302  1463.0  1690.0  2512.0  1410.0  2229
## 6    ssc-let-7f 35432.5 16422 20161.5 22985.5 16512.5 24475.5 19175
```

```r
dim(lma.dfmeanrc)
```

```
## [1] 411 175
```

These dimensions are what would be expected, because there are 411 mature sus scrofa miRNA sequences in miRBase,
and there are 174 animals in the analysis, plus the miRNA column.

Check that the correct miRNA name went with the correct data:


```r
if (sum(names(meanrc)!=lma.dfmeanrc[,1]) != 0) stop ("miRNA names are not the same")

colnames(lma.dfmeanrc)[[1]]<-"miRNA"

if (sum(colnames(lma.dfmeanrc)!=colnames(lmaquant)) != 0) stop ("animal order not the same")

head(lma.dfmeanrc[,1:10])
```

```
##           miRNA     001   002     003     004     005     006   007
## 1    ssc-let-7a 48132.5 23427 28447.5 29860.0 40757.5 29749.5 38799
## 2    ssc-let-7c 32745.0 14987 18144.0 18681.0 34313.0 15294.0 28022
## 3 ssc-let-7d-3p   381.0   192   198.0   269.0   778.0   239.0   774
## 4 ssc-let-7d-5p  4925.0  1938  2511.0  3076.0  3472.0  3276.0  3705
## 5    ssc-let-7e  2811.0  1302  1463.0  1690.0  2512.0  1410.0  2229
## 6    ssc-let-7f 35432.5 16422 20161.5 22985.5 16512.5 24475.5 19175
##       008     009
## 1 31459.5 29857.0
## 2 19512.0 16741.0
## 3   355.0   243.0
## 4  3301.0  3094.0
## 5  1577.0  1523.0
## 6 22094.5 22425.5
```


Set first column of lma.dfmeanrc (miRNA ids) as the row.names:


```r
rownames(lma.dfmeanrc)<-lma.dfmeanrc$miRNA
```

Eliminate column of row names:


```r
lma.dfmeanrc<-lma.dfmeanrc[,-c(1)]
head(lma.dfmeanrc[,1:10])
```

```
##                   001   002     003     004     005     006   007     008
## ssc-let-7a    48132.5 23427 28447.5 29860.0 40757.5 29749.5 38799 31459.5
## ssc-let-7c    32745.0 14987 18144.0 18681.0 34313.0 15294.0 28022 19512.0
## ssc-let-7d-3p   381.0   192   198.0   269.0   778.0   239.0   774   355.0
## ssc-let-7d-5p  4925.0  1938  2511.0  3076.0  3472.0  3276.0  3705  3301.0
## ssc-let-7e     2811.0  1302  1463.0  1690.0  2512.0  1410.0  2229  1577.0
## ssc-let-7f    35432.5 16422 20161.5 22985.5 16512.5 24475.5 19175 22094.5
##                   009   010
## ssc-let-7a    29857.0 21333
## ssc-let-7c    16741.0 12522
## ssc-let-7d-3p   243.0   166
## ssc-let-7d-5p  3094.0  2068
## ssc-let-7e     1523.0   993
## ssc-let-7f    22425.5 16483
```

```r
dim(lma.dfmeanrc)
```

```
## [1] 411 174
```

### 4. Restore the pig IDs in place of the 3-digit codes as the column names of the data frame, for use with the gblup function of gwaR


```r
head(configfile)
```

```
##   pigid code
## 1  1034  001
## 2  1036  002
## 3  1041  003
## 4  1049  004
## 5  1058  005
## 6  1060  006
```

Now I need to substitute the 3-digit code with the pig IDs, ensuring the names stay in the correct order:

Use match function to find positional index and match column names:

The object lma.dfmeanrc has column names that need to be re-named. I have the config file which contains
the current column names and the desired column names. What I am doing in this code is re-ordering the config file 
based on where the config file "code" column matches the position of the lma.dfmeanrc object's column names, then having it return the corresponding value in column "pigid". 

So, when using match, need to have the first argument be the matrix/dataframe you want to change or match, and the second argument be what you want to index it by or match it against. 

"Where does [vector] match in [matrix]?" or "Match the column names of lma.dfmeanrc to the configfile "code" column, then return the corresponding pigid."


```r
configfile[match(colnames(lma.dfmeanrc),configfile$code),"pigid"]
```

```
##   [1] "1034" "1036" "1041" "1049" "1058" "1060" "1080" "1082" "1085" "1091"
##  [11] "1096" "1100" "1107" "1110" "1111" "1113" "1116" "1123" "1134" "1136"
##  [21] "1145" "1147" "1152" "1154" "1158" "1170" "1177" "1179" "1192" "1194"
##  [31] "1197" "1199" "1205" "1207" "1237" "1239" "1240" "1242" "1265" "1267"
##  [41] "1278" "1282" "1291" "1295" "1300" "1304" "1321" "1323" "1423" "1424"
##  [51] "1425" "1426" "1431" "1434" "1435" "1444" "1445" "1449" "1456" "1458"
##  [61] "1482" "1484" "1491" "1493" "1502" "1504" "1510" "1512" "1517" "1523"
##  [71] "1529" "1532" "1533" "1534" "1537" "1543" "1578" "1580" "1589" "1592"
##  [81] "1593" "1594" "1625" "1627" "1638" "1640" "1644" "1646" "1652" "1662"
##  [91] "1669" "1677" "1685" "1687" "1695" "1697" "1746" "1758" "1760" "1776"
## [101] "1778" "1782" "1784" "1785" "1789" "1793" "1798" "1800" "1818" "1819"
## [111] "1820" "1833" "1836" "1839" "1843" "1844" "1879" "1881" "1884" "1886"
## [121] "1889" "1891" "1903" "1904" "1907" "1910" "1914" "1916" "1928" "1930"
## [131] "1965" "1971" "1976" "1980" "1989" "1991" "1999" "2003" "2018" "2020"
## [141] "2022" "2024" "2026" "2027" "2029" "2030" "2064" "2071" "2073" "2076"
## [151] "2094" "2100" "2118" "2119" "2120" "2123" "2131" "2135" "2141" "2143"
## [161] "2152" "2154" "2164" "2168" "2195" "2197" "2229" "2231" "2261" "2263"
## [171] "2297" "2303" "2311" "2317"
```

Assign the column names using match:


```r
colnames(lma.dfmeanrc)<- configfile[match(colnames(lma.dfmeanrc),configfile$code),"pigid"]
head(lma.dfmeanrc[1:10])
```

```
##                  1034  1036    1041    1049    1058    1060  1080    1082
## ssc-let-7a    48132.5 23427 28447.5 29860.0 40757.5 29749.5 38799 31459.5
## ssc-let-7c    32745.0 14987 18144.0 18681.0 34313.0 15294.0 28022 19512.0
## ssc-let-7d-3p   381.0   192   198.0   269.0   778.0   239.0   774   355.0
## ssc-let-7d-5p  4925.0  1938  2511.0  3076.0  3472.0  3276.0  3705  3301.0
## ssc-let-7e     2811.0  1302  1463.0  1690.0  2512.0  1410.0  2229  1577.0
## ssc-let-7f    35432.5 16422 20161.5 22985.5 16512.5 24475.5 19175 22094.5
##                  1085  1091
## ssc-let-7a    29857.0 21333
## ssc-let-7c    16741.0 12522
## ssc-let-7d-3p   243.0   166
## ssc-let-7d-5p  3094.0  2068
## ssc-let-7e     1523.0   993
## ssc-let-7f    22425.5 16483
```

```r
dim(lma.dfmeanrc)
```

```
## [1] 411 174
```

```r
if (sum(colnames(lma.dfmeanrc)!=(configfile$pigid))!=0) stop ("match function did not work correctly")
```

### 5. Round the mean read counts to the nearest integer for use with the voom function


```r
lma.dfmeanrc[1:10,1:6]
```

```
##                   1034   1036     1041     1049     1058     1060
## ssc-let-7a     48132.5  23427  28447.5  29860.0  40757.5  29749.5
## ssc-let-7c     32745.0  14987  18144.0  18681.0  34313.0  15294.0
## ssc-let-7d-3p    381.0    192    198.0    269.0    778.0    239.0
## ssc-let-7d-5p   4925.0   1938   2511.0   3076.0   3472.0   3276.0
## ssc-let-7e      2811.0   1302   1463.0   1690.0   2512.0   1410.0
## ssc-let-7f     35432.5  16422  20161.5  22985.5  16512.5  24475.5
## ssc-let-7g     16700.0   8010   9774.0  10902.0  11420.0  11251.0
## ssc-let-7i     10966.0   4496   5306.0   6397.0   6161.0   5474.0
## ssc-miR-1     307458.0 228189 261470.0 283840.0 109839.0 282455.0
## ssc-miR-100     9539.0   3194   3066.0   3069.0  17474.0   2405.0
```

```r
lma.dfmeanrcround<-round(lma.dfmeanrc)

lma.dfmeanrcround[1:10,1:6]
```

```
##                 1034   1036   1041   1049   1058   1060
## ssc-let-7a     48132  23427  28448  29860  40758  29750
## ssc-let-7c     32745  14987  18144  18681  34313  15294
## ssc-let-7d-3p    381    192    198    269    778    239
## ssc-let-7d-5p   4925   1938   2511   3076   3472   3276
## ssc-let-7e      2811   1302   1463   1690   2512   1410
## ssc-let-7f     35432  16422  20162  22986  16512  24476
## ssc-let-7g     16700   8010   9774  10902  11420  11251
## ssc-let-7i     10966   4496   5306   6397   6161   5474
## ssc-miR-1     307458 228189 261470 283840 109839 282455
## ssc-miR-100     9539   3194   3066   3069  17474   2405
```

The final matrix of rounded, mean read counts needs to have miRNA rownames and Animal ID colnames


```r
head(rownames(lma.dfmeanrcround))
```

```
## [1] "ssc-let-7a"    "ssc-let-7c"    "ssc-let-7d-3p" "ssc-let-7d-5p"
## [5] "ssc-let-7e"    "ssc-let-7f"
```

```r
head(colnames(lma.dfmeanrcround))
```

```
## [1] "1034" "1036" "1041" "1049" "1058" "1060"
```

```r
if (sum(colnames(lma.dfmeanrcround)!= configfile$pigid)!= 0) stop ("rownames do not match pigid")
```

### 6. Extract the 96 libraries from the 174 Bioo Scientific dataset matching the 96 animals with the "lma" selection criteria.


```r
lma96.dfmeanrcround<-lma.dfmeanrcround[,ids]
dim(lma96.dfmeanrcround)
```

```
## [1] 411  96
```

```r
head(lma96.dfmeanrcround)
```

```
##                1034  1036  1041  1049  1058  1060  1080  1082  1085  1091
## ssc-let-7a    48132 23427 28448 29860 40758 29750 38799 31460 29857 21333
## ssc-let-7c    32745 14987 18144 18681 34313 15294 28022 19512 16741 12522
## ssc-let-7d-3p   381   192   198   269   778   239   774   355   243   166
## ssc-let-7d-5p  4925  1938  2511  3076  3472  3276  3705  3301  3094  2068
## ssc-let-7e     2811  1302  1463  1690  2512  1410  2229  1577  1523   993
## ssc-let-7f    35432 16422 20162 22986 16512 24476 19175 22094 22426 16483
##                1096  1100  1107  1110  1111  1113  1116  1123  1134  1136
## ssc-let-7a    35977 21304 43336 46963 41826 41174 36678 45540 59997 34675
## ssc-let-7c    20772 12449 26354 28968 28788 29017 29022 31957 45232 23212
## ssc-let-7d-3p   440   182   389   503   514   449   710   460  1071   411
## ssc-let-7d-5p  3259  1873  5003  5616  4916  4470  3083  5035  6661  3659
## ssc-let-7e     1898   954  3130  3081  2900  2163  2069  2654  3942  2326
## ssc-let-7f    29076 15185 32137 36462 31452 29666 13446 34278 34402 25967
##                1145  1147  1152  1154  1158 1170  1177  1179  1192  1194
## ssc-let-7a    36866 34398 42057 49982 35410 2366 25974 41518 42680 30864
## ssc-let-7c    23731 21229 29962 37824 23684 3077 15389 28463 30503 17962
## ssc-let-7d-3p   356   330   319   750   449  525   339   378   393   256
## ssc-let-7d-5p  3756  3558  4141  5225  3845  175  2610  4251  4195  2588
## ssc-let-7e     2373  2205  2829  2481  2471   75  1345  2686  2324  1716
## ssc-let-7f    28076 27274 29414 29411 26172  576 20826 29466 29266 20800
##                1197  1199  1205  1207  1237  1239  1240  1242  1265  1267
## ssc-let-7a    44394 43819 38598 50950 46225 47411 35912 37383 24522 35894
## ssc-let-7c    31221 31310 25291 37902 32444 32430 24015 25399 14477 25063
## ssc-let-7d-3p   438   407   352   498   569   527   269   326   199   265
## ssc-let-7d-5p  4543  4655  3962  5754  5456  5398  3037  3840  2029  3507
## ssc-let-7e     2746  2749  2181  3820  3060  3045  1752  2440  1034  2036
## ssc-let-7f    31999 32016 28727 36304 34956 34056 27197 25904 20112 24842
##                1278  1282  1291  1295  1300  1304  1321  1323  1423  1424
## ssc-let-7a    15876 42837 36219 30672 28829 28805 24042 30292 25858 31424
## ssc-let-7c    14181 31924 24089 20474 21753 18148 14466 18737 16355 19779
## ssc-let-7d-3p  1033   365   261   268   652   195   209   203   183   236
## ssc-let-7d-5p  1356  4167  3663  3019  2545  2906  2304  2973  2427  3006
## ssc-let-7e      557  2730  2396  1674  1383  1900  1114  1508  1524  1782
## ssc-let-7f    10064 29654 25786 20982 18572 20946 18136 22222 19081 23174
##                1425  1426  1431  1434  1435  1444  1445  1449  1456  1458
## ssc-let-7a    10909 33168 23342 28218 17888 15044 16240 22696 13046 28524
## ssc-let-7c     5973 20260 14482 18147 10798  8526  9964 15666  7201 17704
## ssc-let-7d-3p    84   170   167   212   149    84   101   186    69   147
## ssc-let-7d-5p   826  3037  2314  2058  1582  1233  1421  1378  1004  2518
## ssc-let-7e      349  1968  1254  1413   977   665   728   935   494  1225
## ssc-let-7f     7665 24686 16910 20088 11720  9891 10624 16742  8510 20538
##                1482  1484  1491  1493  1502  1504  1510  1512  1517  1523
## ssc-let-7a    17197 32304 28680 37439 47522 42747 41810 20253 28069 47363
## ssc-let-7c     9342 24491 20807 26015 28812 26738 26126 17415 16928 34465
## ssc-let-7d-3p   116   566   308   379   383   536   323   979   264   506
## ssc-let-7d-5p  1635  2799  2618  3885  4754  4920  4554  1831  2821  4780
## ssc-let-7e      768  1786  1261  2065  3047  2534  2852   996  1365  2816
## ssc-let-7f    12142 12570 21874 27558 35356 32477 31694 10026 21929 35707
##                1529  1532  1533  1534  1537  1543  1578  1580  1589  1592
## ssc-let-7a    29313 26772 35966 27358 42022 44591 48268 35516 30704 37564
## ssc-let-7c    17004 15809 22883 16277 27810 31228 34068 25721 18956 24904
## ssc-let-7d-3p   254   161   310   180   487   507   460   492   252   345
## ssc-let-7d-5p  2946  2586  3729  2315  5000  5502  5521  3304  2918  3832
## ssc-let-7e     1496  1127  2030  1532  2639  2478  3076  1880  1578  2163
## ssc-let-7f    22302 19184 26026 20536 31470 34576 33807 19231 23504 27732
##                1593  1594  1625  1627  1638  1640  1644  1646  1652  1662
## ssc-let-7a    38084 36868 40084 35862 23068 37179 38866 29532 28213 39304
## ssc-let-7c    24512 22762 25993 24232 14980 19671 24475 18429 18530 25909
## ssc-let-7d-3p   369   334   502   292   214   271   306   317   249   417
## ssc-let-7d-5p  4164  3616  4835  3364  2023  3814  3651  2961  2957  3820
## ssc-let-7e     1887  2096  2469  2071  1160  2017  2711  1471  1670  2262
## ssc-let-7f    29918 27982 32046 24424 17376 30385 27878 22126 20542 29883
##                1669  1677  1685  1687  1695  1697
## ssc-let-7a    37340 33902 25400 28632 26621 23326
## ssc-let-7c    23004 20516 14806 15282 13519 13682
## ssc-let-7d-3p   405   265   295   170   293   245
## ssc-let-7d-5p  4027  3215  2537  2731  2630  2136
## ssc-let-7e     2314  1888  1126  1320  1291  1209
## ssc-let-7f    28316 25278 19161 21219 21408 17460
```

```r
colnames(lma96.dfmeanrcround)
```

```
##  [1] "1034" "1036" "1041" "1049" "1058" "1060" "1080" "1082" "1085" "1091"
## [11] "1096" "1100" "1107" "1110" "1111" "1113" "1116" "1123" "1134" "1136"
## [21] "1145" "1147" "1152" "1154" "1158" "1170" "1177" "1179" "1192" "1194"
## [31] "1197" "1199" "1205" "1207" "1237" "1239" "1240" "1242" "1265" "1267"
## [41] "1278" "1282" "1291" "1295" "1300" "1304" "1321" "1323" "1423" "1424"
## [51] "1425" "1426" "1431" "1434" "1435" "1444" "1445" "1449" "1456" "1458"
## [61] "1482" "1484" "1491" "1493" "1502" "1504" "1510" "1512" "1517" "1523"
## [71] "1529" "1532" "1533" "1534" "1537" "1543" "1578" "1580" "1589" "1592"
## [81] "1593" "1594" "1625" "1627" "1638" "1640" "1644" "1646" "1652" "1662"
## [91] "1669" "1677" "1685" "1687" "1695" "1697"
```

```r
sum(colnames(lma96.dfmeanrcround) != ids)
```

```
## [1] 0
```

Check that the data is the same between the subset and the lma.dfmeanrcround dataset:


```r
for (i in colnames(lma96.dfmeanrcround)){
        print(all.equal(lma96.dfmeanrcround[,i], lma.dfmeanrcround[,i]))
}
```

```
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
```

How many miRNAs have 0 expression in this dataset?


```r
sum(rowSums(lma96.dfmeanrcround)==0)
```

```
## [1] 79
```

## Visualize
## Save data
What I am saving here is the rounded, average read counts, and the corresponding covariate data in an .Rdata object


```r
save(lma96.dfmeanrcround, lmapigs, file="../1_lma96_covar_rounded_mean_mature_mirna_expression.Rdata")
```

