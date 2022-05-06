#!/bin/env bash
# author:Stephane Plaisance (VIB-NC), 2019-12-12
# from: https://github.com/BITS-VIB/NGS-Variant-Analysis-training-2020

## Where are we?
workdir="/data/NC_projects/4120_ARuiz_chicken_variants"

## remap reads from SRR17149558 to gallus_gallus GRCg6a.105
cd ${workdir}

# create tmp folder
mkdir -p tmpfiles

# multithreading, adapt to your cpu count
samtoolsthr=24

# how much RAM can java use? - set these to less than your available RAM!
javaopts="-Xms24g -Xmx24g"

#############################################
# PICARD Cleanup & MarkDuplicates
#############################################

reference_fa=/data/biodata/references/galGal6a/Gallus_gallus.GRCg6a.dna.toplevel.fa

infolder=bwa_mappings
outfolder=gatk4_BP

# more records in RAM speeds up when enough RAM is present
recinram=10000000

#################
# function
#################

runall () {

# get bam from the loop
bam=$1
pfx=$(basename ${bam%.bam})

# sort by queryname
cmd="java ${javaopts} -jar $PICARD/picard.jar \
	SortSam \
	I=${bam} \
	O=${outfolder}/${pfx}_qrysrt.bam \
	SO=queryname \
	MAX_RECORDS_IN_RAM=${recinram} \
	TMP_DIR=tmpfiles/"

echo "# ${cmd}"
eval ${cmd}

# mark duplicates
pixdist=100
cmd="java ${javaopts} -jar $PICARD/picard.jar \
	MarkDuplicates \
	I=${outfolder}/${pfx}_qrysrt.bam \
	O=${outfolder}/${pfx}_mrkdup.bam \
	M=${outfolder}/${pfx}_MarkDuplicates.txt \
	ASO=queryname \
	REMOVE_DUPLICATES=false \
	OPTICAL_DUPLICATE_PIXEL_DISTANCE=${pixdist} \
	MAX_RECORDS_IN_RAM=${recinram} \
	TMP_DIR=tmpfiles/"

echo "# ${cmd}"
eval ${cmd}

# sort by coordinate and index
cmd="java ${javaopts} -jar $PICARD/picard.jar \
	SortSam \
	I=${outfolder}/${pfx}_mrkdup.bam \
	O=${outfolder}/${pfx}_mrkdup_srt.bam \
	SO=coordinate \
	CREATE_INDEX=true \
	MAX_RECORDS_IN_RAM=${recinram} \
	TMP_DIR=tmpfiles/"

echo "# ${cmd}"
eval ${cmd}

# fix tags
cmd="java ${javaopts} -jar $PICARD/picard.jar \
	SetNmMdAndUqTags \
	I=${outfolder}/${pfx}_mrkdup_srt.bam \
	O=${outfolder}/${pfx}_mrkdup_srt-tags.bam \
	R=${reference_fa} \
	CREATE_INDEX=true \
	MAX_RECORDS_IN_RAM=${recinram} \
	TMP_DIR=tmpfiles/"

echo "# ${cmd}"
eval ${cmd}

# validate final BAM content
cmd="java ${javaopts} -jar $PICARD/picard.jar \
	ValidateSamFile \
	I=${outfolder}/${pfx}_mrkdup_srt-tags.bam \
	O=${outfolder}/${pfx}_mrkdup_srt-tags.bam_ValidateSamFile.txt \
	R=${reference_fa} \
	M=SUMMARY \
	MO=100 \
	IGNORE_WARNINGS=FALSE \
	MAX_RECORDS_IN_RAM=${recinram} \
	TMP_DIR=tmpfiles/"

echo "# ${cmd}"
eval ${cmd}
}

for bam in ${infolder}/*_rawmappings.bam; do
runall ${bam} &
done

