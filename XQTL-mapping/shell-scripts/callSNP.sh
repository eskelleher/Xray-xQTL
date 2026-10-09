#!/bin/sh
#SBATCH -J callT 
#SBATCH -o callT.o%j
#SBATCH -p batch
#SBATCH -N 1  
#SBATCH -n 16
#SBATCH -t 72:00:00 
#SBATCH --array=1-16
#SBATCH --mail-type=END




module load SAMtools
module load BCFtools



cd /project/kelleher/LLew/Project_xQTL/

#path to the reference genome
ref="/project/kelleher/LLew/Project_xQTL/RIL_data/wgetBAM/R6.founder.bams/dm6.fa"
#make a list of targets for bcftools i.e. chromosome arms
declare -a chrs=("chrX" "chr2L" "chr2R" "chr3L" "chr3R")
mychr=${chrs[$SLURM_ARRAY_TASK_ID - 1]}


# files.txt is a list of bam files to combine
# It should include all the pooled samples AND the DSPR founders as BAMs
#note that file paths willb e different on yoru commuter.
bcftools mpileup -I -d 1000 -t $mychr -a "FORMAT/AD,FORMAT/DP" -f $ref -b RIL_data/wgetBAM/R6.founder.bams/Files/files.txt | bcftools call -mv -Ob > RIL_data/wgetBAM/R6.founder.bams/Xout2/calls.$mychr.bcf  


# the commands below isolate particular fields of the bcftools output and pass them to two separate perl scripts, accuracy.freqtab.pl and accuracy.counttab.pl 
# we did not write these scripts, they were obtained from https://github.com/tdlong/fly_XQTL/tree/main/scripts

# frequencies @ SNPs by sample
bcftools query -e'GT ="./." && QUAL<60' -f'%CHROM %POS %REF %ALT [ %AD{0} %AD{1}] [%GT]\n' RIL_data/wgetBAM/R6.founder.bams/Xout2/calls.$mychr.bcf | grep -v '\.' | perl RIL_data/wgetBAM/R6.founder.bams/Files/accuracy.freqtab.pl > RIL_data/wgetBAM/R6.founder.bams/Xout2/temp.m.$mychr.txt
  
# counts @ SNPs by sample
bcftools query -e 'GT ="./." && QUAL<60' -f'%CHROM %POS %REF %ALT [ %AD{0} %AD{1}] [%GT]\n' RIL_data/wgetBAM/R6.founder.bams/Xout2/calls.$mychr.bcf | grep -v '\.' | perl RIL_data/wgetBAM/R6.founder.bams/Files/accuracy.counttab.pl >RIL_data/wgetBAM/R6.founder.bams/Xout2/temp.countm.$mychr.txt

#change directories to where you freq and count files are simplify the code

cd RIL_data/wgetBAM/R6.founder.bams/Xout2/
# the code below combines the frequency and counts for each SNP from separate files corresponding to separate chromosome ARMS into a single file
#A-D in the file names refer to biological replicates 1-3, C or T indicates control or treated

##first we make the header line
echo -ne "CHROM\tPOS\tfreq_AB8\tfreq_B1\tfreq_B2\tfreq_B3\tfreq_B4\tfreq_B5\tfreq_B6\tfreq_B7\t" > SNP.accuracy.freqXrayA.txt
echo -ne "freq_A.C.100\tfreq_A.T.100\tfreq_B.C.100\tfreq_B.T.100\tfreq_C.C.100\tfreq_C.T.100\tfreq_D.C.100\tfreq_D.T.100\n" >> SNP.accuracy.freqXrayA.txt

##now we append the files for each chromosome
cat temp.chrX.txt >> SNP.accuracy.freqXrayA.txt
cat temp.chr2L.txt >> SNP.accuracy.freqXrayA.txt
cat temp.chr2R.txt >> SNP.accuracy.freqXrayA.txt
cat temp.chr3L.txt >> SNP.accuracy.freqXrayA.txt
cat temp.chr3R.txt >> SNP.accuracy.freqXrayA.txt



echo -ne "CHROM\tPOS\tN_AB8\tN_B1\tN_B2\tN_B3\tN_B4\tN_B5\tN_B6\tN_B7\t" > SNP.accuracy.NXrayA.txt
echo -ne "N_A.C.100\tN_A.T.100\tN_B.C.100\tN_B.T.100\tN_C.C.100\tN_C.T.100\tN_D.C.100\tN_D.T.100\n" >> SNP.accuracy.NXrayA.txt
cat temp.count.chrX.txt >> SNP.accuracy.NXrayA.txt
cat temp.count.chr2L.txt >> SNP.accuracy.NXrayA.txt
cat temp.count.chr2R.txt >> SNP.accuracy.NXrayA.txt
cat temp.count.chr3L.txt >> SNP.accuracy.NXrayA.txt
cat temp.count.chr3R.txt >> SNP.accuracy.NXrayA.txt


# the next step is to do some quality control all on the SNPs to ensure that they are heterozygous in all the founders and can be called in all the founders.