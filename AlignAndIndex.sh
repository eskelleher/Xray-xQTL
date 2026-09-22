#!/bin/bash

#SBATCH -J Index_Dmel 
#SBATCH -o trimBWA_loop_test.o%j
#SBATCH -p batch
#SBATCH -N 2  
#SBATCH -n 32
#SBATCH -t 72:00:00 
#SBATCH --mail-type=END
#SBATCH --mail-user=lgreen7@cougarnet.uh.edu 

#load the tools
module load BWA
module load SAMtools

cd /project/kelleher/LLew/Project_xQTL/UV/data/

ref="/project/kelleher/LLew/Project_xQTL/RIL_data/wgetBAM/R6.founder.bams/dm6.fa" ###location of the bwa index for the reference assembly
dir1="/project/kelleher/LLew/Project_xQTL/UV/data/April/" ### location of the fastq files you want to align
dir2="/project/kelleher/LLew/Project_xQTL/UV/data/April/bam"
shortname="ContX1" ### basename for these samples, used throughout the pipeline

# sample code for 1 fastq file, corresponding to reads from the first control sample
bwa mem -t 8 -M "$ref" "$dir1/SRR38348435_1.fq.gz" "SRR38348435_1.fq.gz" | samtools view -bS - > "$dir2/${shortname}_temp.bam"
samtools sort "$dir2/${shortname}_temp.bam" -o "$dir2/${shortname}.bam"
samtools index "$dir2/${shortname}.bam"
