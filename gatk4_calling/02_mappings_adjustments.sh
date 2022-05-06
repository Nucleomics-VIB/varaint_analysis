#!/bin/env bash
# author:Stephane Plaisance (VIB-NC), 2019-12-12
# from: https://github.com/BITS-VIB/NGS-Variant-Analysis-training-2020

## Where are we?
workdir="/data/NC_projects/4120_ARuiz_chicken_variants/analysis"
cd ${workdir}

# create tmp folder
mkdir -p tmpfiles

# multithreading, adapt to your cpu count
samtoolsthr=24

# how much RAM can java use? - set these to less than your available RAM!
javaopts="-Xms24g -Xmx24g"

# my samtools is here
samtools=$BIOTOOLS/samtools/bin/samtools

#############################################
# PICARD Cleanup & MarkDuplicates
#############################################

outfolder=mappings
mkdir -p ${workdir}/${outfolder}

# reference
reference_fa=reference/Gallus_gallus.GRCg6a.dna.toplevel.fa

# more records in RAM speeds up when enough RAM is present
recinram=10000000

for rawbam in mappings/*_rawmappings.bam; do

# edit in the command below
outpfx=$(basename ${rawbam%_rawmappings.bam})

# sort by queryname
java ${javaopts} -jar $PICARD/picard.jar \
	SortSam \
	I=${outfolder}/${outpfx}_rawmappings.bam \
	O=${outfolder}/${outpfx}_rawmappings_qrysrt.bam \
	SO=queryname \
	MAX_RECORDS_IN_RAM=${recinram} \
	TMP_DIR=tmpfiles/

# mark duplicates
pixdist=2500
java ${javaopts} -jar $PICARD/picard.jar \
	MarkDuplicates \
	I=${outfolder}/${outpfx}_rawmappings_qrysrt.bam \
	O=${outfolder}/${outpfx}_mrkdup.bam \
	M=${outfolder}/${outpfx}_MarkDuplicates.txt \
	ASO=queryname \
	REMOVE_DUPLICATES=false \
	OPTICAL_DUPLICATE_PIXEL_DISTANCE=${pixdist} \
	MAX_RECORDS_IN_RAM=${recinram} \
	TMP_DIR=tmpfiles/

# sort by coordinate and index
java ${javaopts} -jar $PICARD/picard.jar \
	SortSam \
	I=${outfolder}/${outpfx}_mrkdup.bam \
	O=${outfolder}/${outpfx}_mrkdup_srt.bam \
	SO=coordinate \
	CREATE_INDEX=true \
	MAX_RECORDS_IN_RAM=${recinram} \
	TMP_DIR=tmpfiles/

# fix tags
java ${javaopts} -jar $PICARD/picard.jar \
	SetNmMdAndUqTags \
	I=${outfolder}/${outpfx}_mrkdup_srt.bam \
	O=${outfolder}/${outpfx}_mrkdup_srt-tags.bam \
	R=${reference_fa} \
	CREATE_INDEX=true \
	MAX_RECORDS_IN_RAM=${recinram} \
	TMP_DIR=tmpfiles/

# validate final BAM content
java ${javaopts} -jar $PICARD/picard.jar \
	ValidateSamFile \
	I=${outfolder}/${outpfx}_mrkdup_srt-tags.bam \
	O=${outfolder}/${outpfx}_mrkdup_srt-tags.bam_ValidateSamFile.txt \
	R=${reference_fa} \
	M=SUMMARY \
	MO=100 \
	IGNORE_WARNINGS=FALSE \
	MAX_RECORDS_IN_RAM=${recinram} \
	TMP_DIR=tmpfiles/

done
