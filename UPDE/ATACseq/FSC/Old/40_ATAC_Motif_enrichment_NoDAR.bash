#!/bin/bash

#BSUB -L /bin/bash
#BSUB -e /scratch/el/monthly/bdeplanc/pezoldt/Analysis/ATAC_FSC_all/ATAC_FSC_all.err
#BSUB -o /scratch/el/monthly/bdeplanc/pezoldt/Analysis/ATAC_FSC_all/ATAC_FSC_all.out
#BSUB -J Motif[1-4]
#BSUB -M 40000000
#BSUB -R rusage[mem=40000]
#BSUB -n 8
#BSUB -u joern.pezoldt@epfl.ch

# Author by: Joern Pezoldt
# 25.07.2018
# Input:
# a) background bed (e.g. regions not DEG but expressed and not DAR but open)
# Function:
# 1) Find Motifs 
# Note: Run on VitalIT

module use /software/module/;
module add UHTS/Analysis/homer/4.9;

rootdir=/scratch/el/monthly/bdeplanc/pezoldt;
motifdir=$rootdir/Analysis/ATAC_FSC_all/homer/Motif_Input_R/Run_2/TF_reg/BED_Peaks;
backgrounddir=$rootdir/Analysis/ATAC_FSC_all/homer/Motif_Input_R/Run_2/background_genes.txt;
resultdir=$rootdir/Analysis/ATAC_FSC_all/homer/Motif_Input_R/Run_2/TF_reg;
preparseddir=$rootdir/Analysis/ATAC_FSC_all/homer/Motif_Input_R/Run_2;

name=`ls $motifdir | grep '.bed' | sed "${LSB_JOBINDEX}q;d"`;
sample=`ls $motifdir | grep '.bed' | sed "${LSB_JOBINDEX}q;d" | cut -d '.' -f 1`;

echo $sample
echo $name

mkdir -p $resultdir/${sample};

#Find Motifs
findMotifsGenome.pl $motifdir/${name} mm10 $resultdir$/${sample} -size 200 -preparsedDir $preparseddir -bg $backgrounddir
#findMotifs.pl $motifdir/${name} mouse $resultdir/${sample}; -bg $backgrounddir -start -5000 -end 100 -len 8,10,12


