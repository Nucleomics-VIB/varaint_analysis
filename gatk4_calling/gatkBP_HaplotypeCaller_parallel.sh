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

infolder=gatk_BP_recal
outfolder=gatk_BP_HaplotypeCaller
mkdir -p ${outfolder}

# more records in RAM speeds up when enough RAM is present
recinram=10000000

#################
# function
#################

runall () {

# get bam from the loop
bam=$1
pfx=$(basename ${bam%_recal.bam})

#############################################
# GATK4 VARIANT CALLING
#############################################

# parallelize the pair hidden Markov models (pair HMM) process
hmmt=16

java ${javaopts} -jar $GATK/gatk.jar \
	HaplotypeCaller  \
	--input ${bam} \
	--output ${outfolder}/${pfx}.g.vcf.gz \
	--reference ${reference_fa} \
	--emit-ref-confidence GVCF \
	--sample-ploidy 2 \
	--native-pair-hmm-threads ${hmmt} \
	--bam-output ${outfolder}/${pfx}_HC_aligned_reads.bam \
	--tmp-dir tmpfiles/

}

for bam in ${infolder}/*_recal.bam; do
runall ${bam} &
done