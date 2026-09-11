# Xray-xQTL
this respository includes all scripts required for QTL analysis of Xray tolerance in Drosophila, as well as associated RNA-seq data of radiotolerant and radiosensitive strains. The scripts are described more completely below.

# xQTL analysis:

xQTL analysis in the DSPR involves three major stages: 1) SNPs calling and QC,  2) Estimation of haplotype frequencies throughout the genome from pooled sequencing samples based on aforementioned SNP calls, 3) Detection of QTL based on systematic changes in haplotype frequencies in experimental samples when compared to controls at particular genomic positions. The scripts for the of these stages are outlined below. 

# 1) SNP calls and QC

# 2) Haplotype frequency estimates

# 3) QTL mapping and haplotype phasing

LOD.R this script calculates LOD scores throughout the genome based on the following linear model: 

afreq~founder+founder:trt+founder:rep 

where afreq refers to the arcsin transformed haplotype frequency, founder refers to the founder (B1-B8), treatment refers to irradiated (T) or control (C) and replicates is A-D. It also calculates the number of unique haplotypes at each position so that haplotype uncertainty can be easily visualized on the LOD plot.

LOD_haplogroup.R this script offers two important modifications to standard QTL mapping. 1) It compares three window sizes for estimation of haplotype frequencies, 200 Kb (standard), 500 Kb and 1 Mb. It further modifies haplotype frequency data to account for haplotypes that cannot be resolved by combining them into haplogroups. Genotype to phenotype associations are then detected using haplogroup allele frequency changes in experimental vs control populations at particular positions. 

This script also calculates and visualizes that ratio of the haplogroup frequency in replicate experimental and control populations to isolate haplogroups that increase or decrease as a consequence of artificial selection. It investigates how robust the LOD peak is to the removal of particular haplogroups from the data, thereby showing which are most important for driving the LOD peak. 


4) Identification of candidate variants

5) RNAseq
