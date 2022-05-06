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
# GATK4 BAM RECALIBRATION
#############################################

infolder=mappings
outfolder=bam_recalibration
mkdir -p ${outfolder}

dbsnp=reference/ens_dbSNP_150/GRCg6a_dbSNP_150.vcf.gz
reference_fa=reference/Gallus_gallus.GRCg6a.dna.toplevel.fa

# more records in RAM speeds up when enough RAM is present
recinram=10000000

for bamfile in ${infolder}/*_mrkdup_srt-tags.bam; do

outpfx=$(basename ${bamfile%_mrkdup_srt-tags.bam})
recalbamfile=${outpfx}_mrkdup_srt_recal.bam

# compute table before
java ${javaopts} -jar $GATK/gatk.jar \
	BaseRecalibrator \
	-I ${bamfile} \
	-R ${reference_fa} \
	--known-sites ${dbsnp} \
	-O ${outfolder}/${outpfx}_recal_data.table \
	--tmp-dir tmpfiles/

# apply recalibration table
java ${javaopts} -jar $GATK/gatk.jar \
	ApplyBQSR \
	-I ${bamfile} \
	-R ${reference_fa} \
	-bqsr ${outfolder}/${outpfx}_recal_data.table \
	-O ${outfolder}/${recalbamfile} \
	--interval-padding 100 \
	--add-output-sam-program-record \
	--use-original-qualities \
	--static-quantized-quals 10 \
	--static-quantized-quals 20 \
	--static-quantized-quals 30 \
	--tmp-dir tmpfiles/

# compute table after
java ${javaopts} -jar $GATK/gatk.jar \
	BaseRecalibrator \
	-I ${outfolder}/${recalbamfile} \
	-R ${reference_fa} \
	--known-sites ${dbsnp} \
	-O ${outfolder}/${outpfx}_recal_data_after.table \
	--tmp-dir tmpfiles/

# create plots from both tables
java ${javaopts} -jar $GATK/gatk.jar \
	AnalyzeCovariates \
	-before ${outfolder}/${outpfx}_recal_data.table \
	-after ${outfolder}/${outpfx}_recal_data_after.table \
	-plots ${outfolder}/${outpfx}_BQSR_report.pdf \
	-csv ${outfolder}/${outpfx}_BQSR-report.csv \
	--tmp-dir tmpfiles/

# Picard CollectMultipleMetrics on final BAM
java ${javaopts} -jar $PICARD/picard.jar \
	CollectMultipleMetrics \
	I=${outfolder}/${recalbamfile} \
	O=${outfolder}/${outpfx}_multiple_metrics \
	R=${reference_fa} \
	MAX_RECORDS_IN_RAM=${recinram} \
	TMP_DIR=tmpfiles/

done
