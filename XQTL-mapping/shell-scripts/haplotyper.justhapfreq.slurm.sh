#!/bin/bash
#SBATCH --array=1-8
#SBATCH --job-name=haplotyperTrial
#SBATCH --output=haplotyper_%A_%a.out
#SBATCH --error=haplotyper_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your email



module load R

#full path to samples.txt file containing you sample abbreviations and preferred output file names
file="/project/kelleher/LLew/Project_xQTL/RIL_data/wgetBAM/R6.founder.bams/Files/samples.txt"
#full path to the folders where your output files should be written
folder="/project/kelleher/LLew/Project_xQTL/RIL_data/wgetBAM/R6.founder.bams/HapOut"
#full path to file with filtered SNP calls
SNPtable="/project/kelleher/LLew/Project_xQTL/RIL_data/wgetBAM/R6.founder.bams/HapOut/SNP.accuracy.freqXrayF.txt"
#full path to the founder names file
foundernames="/project/kelleher/LLew/Project_xQTL/RIL_data/wgetBAM/R6.founder.bams/Files/ff1.txt"

#reads the codes identifying each pool in samples.txt
pool=$(head -n $SLURM_ARRAY_TASK_ID $file | tail -n 1 | cut -f 1) 
#reads the output files names in samples.txt
OutName=$(head -n $SLURM_ARRAY_TASK_ID $file | tail -n 1 | cut -f 2)

#write the command to your slurm output file
echo "Rscript /project/kelleher/LLew/Project_xQTL/Scripts/haplotyper.justhapfreq2.R $pool $OutName $folder $SNPtable $foundernames"
#calls an R script
Rscript --verbose /project/kelleher/LLew/Project_xQTL/Scripts/haplotyper.justhapfreq2.R $pool $OutName $folder $SNPtable $foundernames
