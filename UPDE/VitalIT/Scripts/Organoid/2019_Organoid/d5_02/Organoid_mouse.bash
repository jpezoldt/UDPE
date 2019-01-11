#!/bin/bash
#BSUB -L /bin/bash
#BSUB -o d502.out                    ##Input here                      ## the standaroutput file
#BSUB -e d502.err                    ##Input here                      ## the standarerror file
#BSUB -J d502                        ##Input here                      ##the job name, change to sample numbers
#BSUB -n 16                                                ## core numbers
#BSUB -R "span[hosts=1]"                                   ## core in the same host, while hosts=2 cross hosts
#BSUB -R "rusage[mem=40000]"                               ## Required memory in MB (at scheduling time)
#BSUB -M 40000000
#BSUB -u jorn.pezoldt@epfl.ch
#BSUB -N

module use /software/module/;
module add UHTS/Aligner/STAR/2.5.3a;
module add UHTS/Analysis/samtools/latest;
module add UHTS/Quality_control/fastqc/0.11.2;
module add Development/java_jdk/latest;
module add UHTS/Analysis/HTSeq/0.6.1;
module add UHTS/Quality_control/RSeQC/2.6.1;
module add R/3.3.2;

#Input here	
name=day_5_c4_S12;																		        
expectedCellNumber=1000;
genomedir=/scratch/el/monthly/bdeplanc/marjan/genomes/mouse;
fastafile=$genomedir/mm_m38r91.fa;
STARIndexdir=$genomedir/STAR_index;
inputGTF=$genomedir/gtf/mm_m38r91.gtf;
#Input here
rootdir=/scratch/el/monthly/bdeplanc/pezoldt/Analysis/Organoid/2019_Organoid/d5_02;				
tooldir=/scratch/el/monthly/bdeplanc/pezoldt/Scripts/Drop-seq_tools-1.13;
dropseqjarfile=$tooldir/jar/dropseq.jar;
#Input here
inputfile=/scratch/el/monthly/bdeplanc/pezoldt/Data/Organoid/2019_Organoid/d5_02/$name;         
outputdir=$rootdir/BAM;
resultdir=$rootdir/Results
QCdir=$resultdir/QC;


echo "Create directories"
mkdir -p $resultdir/Coverage
mkdir -p $outputdir/
mkdir -p $QCdir/

echo "Processing Sample" $name " following DropSeq pipeline...";
echo "1. Fastq to SAM using Picard"
java -jar $tooldir/3rdParty/picard/picard.jar FastqToSam F1=${inputfile}_R1_001.fastq.gz F2=${inputfile}_R2_001.fastq.gz QUALITY_FORMAT=Standard O=$outputdir/$name.bam SM=$name

echo "2. Filtering"
echo "2.1. Tag Cells (12bp Barcode)"
java -jar $dropseqjarfile TagBamWithReadSequenceExtended SUMMARY=$resultdir/cellular_tag_summary.txt BASE_RANGE=1-12 BASE_QUALITY=10 BARCODED_READ=1 DISCARD_READ=false TAG_NAME=XC NUM_BASES_BELOW_QUALITY=1 INPUT=$outputdir/$name.bam OUTPUT=$outputdir/$name.cTag.bam
echo "2.2. Tag Molecules (8bp UMI)"
java -jar $dropseqjarfile TagBamWithReadSequenceExtended SUMMARY=$resultdir/molecular_tag_summary.txt BASE_RANGE=13-20 BASE_QUALITY=10 BARCODED_READ=1 DISCARD_READ=true TAG_NAME=XM NUM_BASES_BELOW_QUALITY=1 INPUT=$outputdir/$name.cTag.bam OUTPUT=$outputdir/$name.cTag.mTag.bam
echo "2.3. Filter Bad quality reads"
java -jar $dropseqjarfile FilterBAM TAG_REJECT=XQ INPUT=$outputdir/$name.cTag.mTag.bam OUTPUT=$outputdir/$name.cTag.mTag.filtered.bam
echo "2.4. Trim SMART adapter in 5' sequence of reads (at least 5bp)"
java -jar $dropseqjarfile TrimStartingSequence OUTPUT_SUMMARY=$resultdir/adapter_trimming_summary.txt SEQUENCE=AAGCAGTGGTATCAACGCAGAGTGAATGGG MISMATCHES=0 NUM_BASES=5 INPUT=$outputdir/$name.cTag.mTag.filtered.bam OUTPUT=$outputdir/$name.cTag.mTag.filtered.adapTrimmed.bam
echo "2.5. Trim PolyA sequences in 3' sequence of reads (at least 6bp)"
java -jar $dropseqjarfile PolyATrimmer OUTPUT_SUMMARY=$resultdir/polyA_trimming_summary.txt MISMATCHES=0 NUM_BASES=6 INPUT=$outputdir/$name.cTag.mTag.filtered.adapTrimmed.bam OUTPUT=$outputdir/$name.cTag.mTag.filtered.adapTrimmed.polyATrimmed.bam
echo "3. Revert back from BAM to Fastq file"
java -jar $tooldir/3rdParty/picard/picard.jar SamToFastq INPUT=$outputdir/$name.cTag.mTag.filtered.adapTrimmed.polyATrimmed.bam FASTQ=$inputfile.fastq.gz
echo "4. Alignment"
echo "4.1. Mapping with STAR (unsorted)"
STAR --runMode alignReads --runThreadN 4 --genomeDir $STARIndexdir --outFilterMultimapNmax 1 --readFilesCommand zcat --outSAMtype BAM Unsorted --outFileNamePrefix $outputdir/ --readFilesIn $inputfile.fastq.gz
echo "4.2. Sorting by read name (required)"
java -jar $tooldir/3rdParty/picard/picard.jar SortSam INPUT=$outputdir/Aligned.out.bam OUTPUT=$outputdir/$name.sortedByName.bam SO=queryname
echo "4.3. Add tags from unaligned BAM to the STAR-aligned BAM file"
java -jar $tooldir/3rdParty/picard/picard.jar MergeBamAlignment UNMAPPED_BAM=$outputdir/$name.cTag.mTag.filtered.adapTrimmed.polyATrimmed.bam ALIGNED_BAM=$outputdir/$name.sortedByName.bam REFERENCE_SEQUENCE=$fastafile PAIRED_RUN=false OUTPUT=$outputdir/$name.merged.bam
echo "4.4. Add GE (Gene Exon) tags from the GTF file to the BAM file when read overlaps a gene exon"
java -jar $dropseqjarfile TagReadWithGeneExon INPUT=$outputdir/$name.merged.bam OUTPUT=$outputdir/$name.merged.GETagged.bam ANNOTATIONS_FILE=$inputGTF TAG=GE
echo "4.5. Detect and Correct Bead Synthesis Errors (NUM_BARCODES parameters should be roughly 2x expected number of cells)"
java -jar $dropseqjarfile DetectBeadSynthesisErrors INPUT=$outputdir/$name.merged.GETagged.bam O=$outputdir/$name.merged.GETagged.cleaned.bam OUTPUT_STATS=$resultdir/synthesis_stats.txt CREATE_INDEX=true SUMMARY=$resultdir/synthesis_stats_summary.txt NUM_BARCODES=$(($expectedCellNumber * 2)) PRIMER_SEQUENCE=AAGCAGTGGTATCAACGCAGAGTAC
echo "5. Digital Gene Expression and final stats/cleaning"
echo "5.1. Generate Histogram of reads per cell barcodes (to select relevant cell number cutoff)"
java -jar $dropseqjarfile BAMTagHistogram INPUT=$outputdir/$name.merged.GETagged.cleaned.bam OUTPUT=$resultdir/histogram_reads_cell_barcodes.txt.gz TAG=XC
#Rscript $tooldir/Histogram.DropSeq.R $resultdir/histogram_reads_cell_barcodes $(($expectedCellNumber * 4))
echo "5.2. Adapt 'expectedCellNumber' variable if needed"
echo "Current value =" $expectedCellNumber
echo "5.2. Generate UMI Matrix (use at least one of MIN_NUM_GENES_PER_CELL, MIN_NUM_TRANSCRIPTS_PER_CELL, NUM_CORE_BARCODES or CELL_BC_FILE)"
java -jar $dropseqjarfile DigitalExpression INPUT=$outputdir/$name.merged.GETagged.cleaned.bam OUTPUT=$resultdir/UMI.counts.dge.txt.gz SUMMARY=$resultdir/dge_summary.txt NUM_CORE_BARCODES=$expectedCellNumber
echo "5.3. Clean intermediate files"
mv $outputdir/Log.final.out $resultdir/alignment_summary.txt
find $outputdir/ -type f -not -name '*.cleaned.ba*' -exec rm {} \;
mv $outputdir/$name.merged.GETagged.cleaned.bam $outputdir/$name.bam
mv $outputdir/$name.merged.GETagged.cleaned.bai $outputdir/$name.bam.bai
echo "5.4. Generate fastQC for fastq file and final BAM"
fastqc --outdir=$QCdir/ $inputfile.fastq.gz
fastqc --outdir=$QCdir/ $outputdir/$name.bam
#echo "5.5. Generate Coverage file"
#geneBody_coverage.py -r $bedfile -i $outputdir/$name.bam -o $resultdir/Coverage
