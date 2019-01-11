#!/bin/bash

#BSUB -L /bin/bash
#BSUB -e /scratch/el/monthly/bdeplanc/pezoldt/Analysis/ATAC_FSC_all/std_err_out/ATAC_FSC_all_v46_sigi_bgPeaks_Ext.err
#BSUB -o /scratch/el/monthly/bdeplanc/pezoldt/Analysis/ATAC_FSC_all/std_err_out/ATAC_FSC_all_v46_sigi_tgPeaks_Ext.out
#BSUB -J Motif[1-8]
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
#homer 4.2 is not on vital-it anymore
module add UHTS/Analysis/homer/4.9;


parameter=sigi_v46_bgPeaks_Ext_;
rootdir=/scratch/el/monthly/bdeplanc/pezoldt;
motifdir=$rootdir/Analysis/ATAC_FSC_all/homer/Motif_Input_R/Run_4_Ext50bp/BED_peaks;
backgrounddir=$rootdir/Analysis/ATAC_FSC_all/homer/Motif_Input_R/Run_4_Ext50bp/background_all_peaks_minus_input_Extension50bp.bed;
resultdir=$rootdir/Analysis/ATAC_FSC_all/homer/Motif_Input_R/Run_4_Ext50bp/Output_Regions;
genomedir=$rootdir/Data/Genome_mouse/GRCm38.87/Mus_musculus.GRCm38.87.dna.primary_assembly.fa
preparseddir=$rootdir/tmp;
name=`ls $motifdir | grep '.bed' | sed "${LSB_JOBINDEX}q;d"`;
sample=`ls $motifdir | grep '.bed' | sed "${LSB_JOBINDEX}q;d" | cut -d '.' -f 1`;

echo $sample
echo $name

#mkdir -p $resultdir/${sample};

#Find Motifs
#findMotifs.pl $motifdir/${name} mouse $resultdir/${sample} -bg $backgrounddir -start -2000 -end 200 -len 8,10,12
#findMotifsGenome.pl $motifdir/${name} mm10 $resultdir/${sample} -preparse -preparsedDir $preparseddir
findMotifsGenome.pl $motifdir/${name} $genomedir $resultdir/$parameter${sample} -bg $backgrounddir -size given -len 10
#$parameter
#$genomedir


