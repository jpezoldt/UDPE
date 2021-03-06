#!/bin/bash

#BSUB -L /bin/bash
#BSUB -e /scratch/el/monthly/bdeplanc/pezoldt/Analysis/WGBS_FSC/std_err_out/WGBS_sigi_bgHOMER_mask_40cpg_split_thresh.err
#BSUB -o /scratch/el/monthly/bdeplanc/pezoldt/Analysis/WGBS_FSC/std_err_out/WGBS_sigi_bgHOMER_mask_40cpg_split_thresh.out
#BSUB -J Motif[1-16]
#BSUB -M 40000000
#BSUB -R rusage[mem=40000]
#BSUB -n 16
#BSUB -u joern.pezoldt@epfl.ch
# Author by: Joern Pezoldt
# 25.09.2018
# Input:
# a) bed with TSS regions of genes associated with open chromatin
# b) bed with DMRs
# Function:
# 1) Find Motifs 
# Note: Run on VitalIT

module use /software/module/;
#homer 4.2 is not on vital-it anymore
module add UHTS/Analysis/homer/4.6;

parameter=sigi_bgHOMER_mask_40cpg_split_thresh_;
rootdir=/scratch/el/monthly/bdeplanc/pezoldt;
motifdir=$rootdir/Analysis/WGBS_FSC/Regions_of_Interest/Run3_CpGRegions_40;
backgrounddir=$rootdir/Analysis/ATAC_FSC_all/homer/Motif_Input_R/Run_6_noExt/background_all_TSS_minus_input_genes.bed;
resultdir=$rootdir/Analysis/WGBS_FSC/Homer_Output/Run3_CpGRegions_40;
genomedir=$rootdir/Data/Genome_mouse/GRCm38.87/Mus_musculus.GRCm38.87.dna.primary_assembly.fa;
preparseddir=$rootdir/tmp;

name=`ls $motifdir | grep '.bed' | sed "${LSB_JOBINDEX}q;d"`;
sample=`ls $motifdir | grep '.bed' | sed "${LSB_JOBINDEX}q;d" | cut -d '.' -f 1`;

echo $sample
echo $name

#mkdir -p $resultdir/${sample};

#Find Motifs
#findMotifs.pl $motifdir/${name} mouse $resultdir/${sample} -bg $backgrounddir -start -2000 -end 200 -len 8,10,12
#findMotifsGenome.pl $motifdir/${name} mm10 $resultdir/${sample} -preparse -preparsedDir $preparseddir
findMotifsGenome.pl $motifdir/${name} $genomedir $resultdir/$parameter${sample} -size given -len 10
# -bg $backgrounddir 
# -mask -cpg

