#!/bin/bash

# Author by: Joern Pezoldt
# 23.08.2018
# Inputs:
# a) PWM file obtained via HOMER
# b) FASTA file for the regions piped into HOMER that allowed the identification of motifs and their respective PWM (see a)
# c) 
# Function:
# 1) Make PWM file fit to gTEMPLE
# 2) Transform BED file into FASTA according to gTEMPLE format
# 3) Identfiy which peak contains which motif and thereby allow to track which TF influences which gene 
# Run on FAMEAUX

#Set Name of Analysis, will be used for:
# a) Task list gTEMPLE
# b) Accessing Motif and BED file
JOB_ID="600bp_Closed_NoDEG_SPF"
#!!!Note: Region files still need to be changed from chr1 to 1, example folder given
JOB_BED="Closed_NoDEG_SPF_noCHR.bed"
#Names output files
PWM_FILE="gTEMPLE_PWM_$JOB_ID"
REGIONS_FILE="gTEMPLE_REGIONS_$JOB_ID"
SEQUENCE_FILE="gTEMPLE_SEQUENCE_$JOB_ID"

#Set Directories
rootdir=/home/pezoldt/NAS2/pezoldt;
fileFASTAGENOME=/data/genome/mus_musculus/GRCm38.91/Mus_musculus.GRCm38.91.dna.primary_assembly.fa;
#set dir for motif file
dirPWM=$rootdir/Analysis/ATACseq/ATAC_FSC_all/Motif/Homer_Output/allHomerResults/$JOB_ID/homerMotifs.motifs10;
dirBED=$rootdir/Analysis/ATACseq/ATAC_FSC_all/Motif/Homer_Output/gTEMPLE/Region_file_Test/$JOB_BED;
dirRESULT=$rootdir/Analysis/ATACseq/ATAC_FSC_all/Motif/Homer_Output/gTEMPLE;

# 1) Make PWM file fit to gTEMPLE 
PWM_file=`cat $dirPWM/homerMotifs.motifs10`
while IFS=$'>' read -r -a myArray
	do
	LEN=$(echo ${#myArray[0]})
	if [ $LEN -lt 10 ]; then
		IFS=$'\t' read -r -a content <<< "${myArray[1]}"
		MOTIF="${content[0]}"
		START=">"
		PVALUE="pval=0.001"
		TAB="\t"
		NEW_HEADER="$START$MOTIF$TAB$PVALUE"
		echo -e "$NEW_HEADER"
	else
		echo "${myArray[0]}"
	fi
done < $dirPWM > $dirRESULT/$PWM_FILE

# 2) Transform BED file into FASTA according to gTEMPLE sequence format
seq_gT=`bedtools getfasta -name+ -fi $fileFASTAGENOME -bed $dirBED`
for line in $seq_gT;
	do
	IFS=: read -r -a content <<< "$line"
	#Determine the length of the first entry in content array
	# if length is shorter 30 it is a header line
	# if length is longer than 30 it is a sequence
	LEN=$(echo ${#content[0]})
	#echo $LEN
	if [ $LEN -lt 30 ]; then
        #echo ${content[2]}
		CHR=`cut -d ':' -f3 <<< $line`
		START=`cut -d ':' -f4 <<< $line | cut -d '-' -f1`
		END=`cut -d ':' -f4 <<< $line | cut -d '-' -f2`
		ID=`cut -d ':' -f1 <<< $line | cut -d '>' -f2`
		input=">chr"
		underscore="_"
		#echo $ID
		#echo $CHR
		#echo $START
		#echo $END
		NEW_HEADER="$input$CHR$underscore$START$underscore$END$underscore$ID"
		echo $NEW_HEADER
	else
        echo $line
	fi
done > $dirRESULT/$SEQUENCE_FILE

#Add the topline header
double_slash="//"
top_line="$double_slash$JOB_ID"
echo "$JOB_ID" | sed 's/ /\\ /g'
sed -i $'1 i\\\n'$top_line'' $dirRESULT/$SEQUENCE_FILE

#Generate Task/Region File for gTEMPLE
for line in $seq_gT;
	do
	IFS=: read -r -a content <<< "$line"
	#Determine the length of the first entry in content array
	# if length is shorter 30 it is a header line
	# if length is longer than 30 it is a sequence
	LEN=$(echo ${#content[0]})
	#echo $LEN
	if [ $LEN -lt 30 ]; then
        #echo ${content[2]}
		COMPARISON=$JOB_ID
		#echo $JOB_ID
		CHR=`cut -d ':' -f3 <<< $line`
		#echo $CHR
		START=`cut -d ':' -f4 <<< $line | cut -d '-' -f1`
		#echo $START
		END=`cut -d ':' -f4 <<< $line | cut -d '-' -f2`
		#echo $END
		ID=`cut -d ':' -f1 <<< $line | cut -d '>' -f2`
		#echo $ID
		TAB="\t"
		PWMS="all"
		CHR_both="chr$CHR"
		#echo $ID
		#echo $CHR
		#echo $START
		#echo $END
		TASK_i="$COMPARISON$TAB$CHR_both$TAB$START$TAB$END$TAB$ID$TAB$PWMS"
		echo -e "$TASK_i"
	fi
done > $dirRESULT/$REGIONS_FILE

#Add headers to Region files
sed -i $'1 i\\\nset\tchr\tstart\tstop\tid\tpwms' $dirRESULT/$REGIONS_FILE

#Run Temple
cd /home/pezoldt/NAS2/pezoldt/Software/Java/gTEMPLE
java -jar gTemple_v1.0.jar $dirRESULT/$PWM_FILE $dirRESULT/$SEQUENCE_FILE $dirRESULT/$REGIONS_FILE 2 0.5 1
