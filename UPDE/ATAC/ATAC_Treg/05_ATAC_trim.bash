#!/bin/bash

#BSUB -L /bin/bash
#BSUB -e /scratch/el/monthly/bdeplanc/pezoldt/Analysis/ATAC_denovo_Treg/ATAC_denovo_Treg.err
#BSUB -o /scratch/el/monthly/bdeplanc/pezoldt/Analysis/ATAC_denovo_Treg/ATAC_denovo_Treg.out
#BSUB -J Treg[1-16]
#BSUB -M 40000000
#BSUB -R rusage[mem=40000]
#BSUB -n 16
#BSUB -u joern.pezoldt@epfl.ch

# Author: Joern Pezoldt
# 17.07.2018
# Function:
# 1) Trims adapter sequences of ATAC-seq
# Note: Run on VitalIT

module use /software/module/;
module add UHTS/Quality_control/cutadapt/1.8;

#Manually change directories
rootdir=/scratch/el/monthly/bdeplanc/pezoldt;
fastqdir=$rootdir/Data/ATAC_denovo_Treg;
fastqdir_trim=$rootdir/Data/ATAC_denovo_Treg/Trimmed;

#Trimming folder
name=`ls $fastqdir/R1 | grep '.fastq.gz$' | sed "${LSB_JOBINDEX}q;d" | cut -d '_' -f 1-5`;
sample=`echo $name |cut -d '_' -f 1-2`;

echo $name
echo $sample


#make directories for trimmed results
mkdir -p $fastqdir/Trimmed;
mkdir -p $fastqdir_trim/R1;
mkdir -p $fastqdir_trim/R2;
fastqdir_trim=$rootdir/Data/ATAC_denovo_Treg/Trimmed;

#trims low quality bases with <30 previous to adapter trimming
#keep trimmed reads with at least length 30
cutadapt -q 30,30 -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -A CTGTCTCTTATACACATCTGACGCTGCCGACGA -m 30 -o $fastqdir_trim/R1/${sample}_1.fastq.gz -p $fastqdir_trim/R2/${sample}_2.fastq.gz $fastqdir/R1/${sample}_1.fastq.gz $fastqdir/R2/${sample}_2.fastq.gz
cutadapt -q 30,30 -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -A CTGTCTCTTATACACATCTGACGCTGCCGACGA -m 30 -o $fastqdir_trim/R1/ATAC_pLNSPF4_1.fastq.gz -p $fastqdir_trim/R2/ATAC_pLNSPF4_2.fastq.gz $fastqdir/R1/ATAC_pLNSPF4_1.fastq.gz $fastqdir/R2/ATAC_pLNSPF4_2.fastq.gz

