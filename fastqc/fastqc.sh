#!/bin/bash 
#SBATCH -A nereus                ### Account
#SBATCH --partition=compute         ### Partition
#SBATCH --job-name=fastqc         ### Job Name
#SBATCH --time=5:00:00        ### WallTime
#SBATCH --nodes=1              ### Number of Nodes
#SBATCH --ntasks=1             ### Number of tasks per array job
##SBATCH --array=0-3           ### Array index

#SBATCH --mail-type=END              ### Mail events (NONE, BEGIN, END, FA$
#SBATCH --mail-user=ssnyder3@uoregon.edu  ### Where to send mail
#SBATCH --cpus-per-task=8            ### Number of CPU cores per task
#SBATCH --requeue


##        DESCRIPTION      ##
# fastqc for quality control on raw sequence data.
# multiqc to combine individual fastqc files into a comprehensive report

##     TO USE    ## 
# add path to raw fastq files

# load modules - FASTQC 
# note: loading modules for both fastqc and multiqc  at once creates problems on Talapas 
module load fastqc/0.11.5

mkdir fastqc_out
mkdir multiqc_out
OUTDIR=fastqc_out

#input fastq dir 
FASTQDIR=/home/ssnyder3/nereus/aging_rnaseq/all_raw_fastqs

# generate fastqc files
#/usr/bin/time fastqc -o $OUTDIR -t 6 $FASTQDIR/*.fastq.gz


# load modules for multiQC
module load eb-hide/1  icc/2017.1.132-GCC-6.3.0-2.27  impi/2017.1.132
module load MultiQC/1.3-Python-3.6.1

# generate multiQC report to summarize 
/usr/bin/time multiqc $OUTDIR -o multiqc_out

