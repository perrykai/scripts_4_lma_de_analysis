#' **Script:** `4_download_targetscan_script_databases.sh`
#' 
#' **Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/4_lma_de_analysis/scripts`
#' 
#' **Date:**  11/28/16
#' 
#' **Input File Directory:**  `http://www.targetscan.org/vert_71/vert_71_data_download/targetscan_70.zip`
#' 
#' **Input File(s):** `targetscan_70.zip`
#' 
#' **Output File Directory:** 
#' 
#' 1. `/mnt/home/perrykai/targetscan_70`
#' 2. `/mnt/research/pigeqtl/analyses/microRNA/reference_sequences`
#' 
#' **Output File(s):** 
#' 
#' 1. `miR_Family_info_sample.txt`
#' 2. `README_70.txt`
#' 3. `targetscan_70_output.txt`
#' 4. `targetscan_70.pl`
#' 5. `targetscan_70.zip`
#' 6. `UTR_Sequences_sample.txt`
#' 7. ``
#' 8. ``
#' **Table of contents:**
#'
#' 1. [Objectives](#objectives)
#' 2. [Analysis](#analysis)
#' 
#' ## Objectives
#' 
#' The objective of this script is to download the targetscan perl script and the necessary databases for use in target prediction 
#' of the significant miRNAs from the lma DE analysis. 
#' Note: On 12/6/16 I downloaded the Gene_info.txt file to match the significant DE miRNA targets' transcript IDs to their gene symbols

#' ## Analysis
cd /mnt/home/perrykai/

wget http://www.targetscan.org/vert_71/vert_71_data_download/targetscan_70.zip

# --2016-11-28 15:05:17--  http://www.targetscan.org/vert_71/vert_71_data_download/targetscan_70.zip
# Resolving www.targetscan.org... 18.4.1.128
# Connecting to www.targetscan.org|18.4.1.128|:80... connected.
# HTTP request sent, awaiting response... 200 OK
# Length: 35766 (35K) [application/zip]
# Saving to: “targetscan_70.zip”

# 100%[===============================================================================>] 35,766      --.-K/s   in 0.09s

# 2016-11-28 15:05:18 (368 KB/s) - “targetscan_70.zip” saved [35766/35766]

mkdir targetscan_70

mv targetscan_70.zip targetscan_70

unzip targetscan_70.zip

ls -ltrh

# total 150K
# -rw-r--r-- 1 perrykai gene-expression  90K Dec 10  2008 UTR_Sequences_sample.txt
# -rw-r--r-- 1 perrykai gene-expression  405 Aug 14  2015 miR_Family_info_sample.txt
# -rw-r--r-- 1 perrykai gene-expression  17K Aug 14  2015 targetscan_70_output.txt
# -rwxr-xr-x 1 perrykai gene-expression  30K Sep 28  2015 targetscan_70.pl
# -rw-r--r-- 1 perrykai gene-expression 4.0K Sep 28  2015 README_70.txt
# -rw-rw-r-- 1 perrykai gene-expression  35K Sep 28  2015 targetscan_70.zip


cat README_70.txt

# INTRODUCTION

# The Perl script targetscan_70.pl identifies miRNA targets and then determines whether a given target is conserved or not across a given set of species.
# The TargetScan 7.0 prediction code produces essentially the same output as the previous version (targetscan_60.pl)
# except in the way that a group of heterogeneous aligned sites (combinations of at least two of the three site types) is classified,
# TargetScan 7.0 defines a "group type" in these cases to be a combination of site types, which allows subgrouping by species that share the same site type.

# The script takes two input files
# 	1) A tab-delimited file that lists the miRNA seed sequences and the species in which they are present.
# 	2) A tab-delimited multiple sequence alignment of the 3' UTRs of genes from the desired species.

# In this directory we are providing samples of both the above files and the script uses them by default.

# The sample files are: UTR_sequences_sample.txt and miR_Family_info_sample.txt
# If you wish to generate these files from the complete data available for download,
# 	run the commands shown below to convert them to the correct format:


# FILE FORMATS

# The format of the input files is important for the script to work correctly.

# Each line of the miRNA seed sequence file consists of 3 tab separated entries
# 1) Name of the miRNA family
# 2) The 7 nucleotide long seed region sequence.
# 3) List (semicolon-delimited) of Species IDs of this miRNA family (which should match species IDs in UTR input file)

# Each line of the alignment file consists of 3 tab separated entries
# 1) Gene symbol or transcript ID
# 2) Species/taxonomy ID (which should match species IDs in miRNA input file)
# 3) Sequence: all uppercase except for regions that overlap an ORF (which will be flagged as such and ignored for subsequent PCT calculation)

# If you wish to generate this files from the complete data available by download,
# 	run these commands to convert it to the correct format (and without a header):
# sed '1,1d' UTR_Sequences.txt | cut -f1,4,5 > UTR_sequences_all.txt


# EXECUTION

# The script can be executed in 3 different ways:
# 1) Running the script without any arguments (./targetscan_70.pl) will print out a help screen.
# 1) Running the script without the '-h' flag (./targetscan_70.pl -h) will print out a formats of input files.
# 2) Running the script with input filenames and output file will perform the analysis. Ex:
# 	./targetscan_70.pl miR_Family_info_sample.txt UTR_Sequences_sample.txt targetscan_70_output.txt

# OUTPUT FILES

# In this folder is a sample output file called "targetscan_70_output.txt".
# The output file contains several tab separated entries per line:

# 	The sample output file has a headers that names each column
# 	GeneID - name/ID of gene (from UTR input file)
# 	miRNA_family_ID - name/ID of miRNA family (from miRNA input file)
# 	species_ID - name/ID of species (from UTR input file)
# 	MSA_start - starting position of site in aligned UTR (counting gaps)
# 	MSA_end - ending position of site in aligned UTR (counting gaps)
# 	UTR_start - starting position of site in UTR (not counting gaps)
# 	UTR_end - ending position of site in UTR (not counting gaps)
# 	Group_ID - ID (number) of site(s) (same gene, same miRNA) that overlap
# 	Site_type - type of site in this species (m8 [7mer-m8], 1a [7mer-1A], or m8:1a [8mer])
# 	miRNA in this species - if "x", then this miRNA has been annotated in this species
# 	Group_type - type of this group of sites; if 'Site_type' in a 'Group_ID' is heterogeneous, "weakest" type of the group is used
# 	Species_in_this_group - list of species names/IDs in which this site is found
# 	Species_in_this_group_with_this_site_type - for hetergeneous groups only
# 	ORF_overlap - If site in the UTR sequence is lowercase (indicating ORF overlap), this will be set to 1.  Typical UTR sites have a value of 0.

# NOTES

# This script was designed on a Linux platform. While running this script on Windows or Mac platforms, make sure to call the native perl binary.


# QUESTIONS/SUGGESTIONS:

# Please direct all correspondence to wibr-bioinformatics@wi.mit.edu


#' 
#' Now to download the miRNA dataset and the mRNA dataset from the targetscan website:
#' 
cd /mnt/research/pigeqtl/analyses/microRNA/reference_sequences/
#' First, the miRNA info:
wget http://www.targetscan.org/vert_71/vert_71_data_download/miR_Family_Info.txt.zip 

#' Second, the mRNA 3'UTR info:
wget http://www.targetscan.org/vert_71/vert_71_data_download/UTR_Sequences.txt.zip

mkdir targetscan_databases

mv miR_Family_Info.txt.zip targetscan_databases
mv UTR_Sequences.txt.zip targetscan_databases

cd /mnt/research/pigeqtl/analyses/microRNA/reference_sequences/targetscan_databases

#' Third, the mRNA "gene info" file:
wget http://www.targetscan.org/vert_71/vert_71_data_download/Gene_info.txt.zip

#' Unzip the full 3'UTR dataset (6Gb) to format it correctly and extract the human information:
unzip UTR_Sequences.txt.zip

#' Format the file correctly for use in targetscan:
sed '1,1d' UTR_Sequences.txt | cut -f1,4,5 > UTR_sequences_all.txt

#' Remove the giant unzipped file:
rm UTR_Sequences.txt

#' The species ID (column 2) for homo sapiens is 9606; use this to extract only the human sequences
awk '$2 == 9606' UTR_sequences_all.txt > UTR_sequences_hsa.txt

#' Then remove the giant files that are the 3' UTR information (zipped file and UTR_sequences_all.txt)
rm UTR_Sequences.txt.zip
rm UTR_sequences_all.txt

#' Unzip the miRNA information:
unzip miR_Family_Info.txt

head miR_Family_Info.txt
# miR family	Seed+m8	Species ID	MiRBase ID	Mature sequence	Family Conservation?	MiRBase Accession
# miR-7119-3p	AAAAACA	10090	mmu-miR-7119-3p	AAAAAACACCGUUUCCUCCAG	-1	MIMAT0028136
# miR-548bd-3p	AAAAACC	9544	mml-miR-548b	CAAAAACCUCAAUUGCUUUUGU	-1	MIMAT0006411
# miR-548bd-3p	AAAAACC	9544	mml-miR-548d-3p	CAAAAACCACAAUUUCUUUUGC	-1	MIMAT0006414
# miR-2285u	AAAAACC	9913	bta-miR-2285u	GAAAAACCCGAACGAACUUU	-1	MIMAT0031094
# miR-548h-3p	AAAAACU	9544	mml-miR-548h-3p	CAAAAACUGCAGUGACUUCUGU	-1	MIMAT0028257
# miR-2284aa/2284z	AAAAAGU	9913	bta-miR-2284aa	AAAAAAGUUUGUUUGGGUUUU	-1	MIMAT0025560
# miR-2284aa/2284z	AAAAAGU	9913	bta-miR-2284z	AAAAAAGUUUGUUUGGGUUUUU	-1	MIMAT0025559
# miR-548c	AAAAAUC	9598	ptr-miR-548c	CAAAAAUCUCAAUUACUUUUGC	-1	MIMAT0008220
# miR-2285x	AAAAAUC	9913	bta-miR-2285x	GAAAAAUCUGAAUGAACUUUUGG	-1	MIMAT0029946

#' Now extract rows in which the third column contains the human ID (9606), and retain only the first three columns.
#'
#' This is all the data that is needed for targetscan: miRNA family, seed sequence, species ID 
awk '$3 == 9606' miR_Family_Info.txt | cut -f1-3 > miR_Family_Info_hsa.txt

head miR_Family_Info_hsa.txt

# miR-9500	AGGGAAG	9606
# miR-8485	ACACACA	9606
# miR-8088	CUCGGUA	9606
# miR-8087	AAGACUU	9606
# miR-8086	GCUAGUC	9606
# miR-8084	AAUACUA	9606
# miR-8083	AGGACUU	9606
# miR-8081	UUGAGUC	9606
# miR-8080	AAGGACA	9606
# miR-8078	GUCUAGG	9606

unzip Gene_info.txt.zip
rm Gene_info.txt.zip