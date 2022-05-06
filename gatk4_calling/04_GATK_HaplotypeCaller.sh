#!/bin/env bash
# author:Stephane Plaisance (VIB-NC), 2019-12-12
# from: https://github.com/BITS-VIB/NGS-Variant-Analysis-training-2020

## Where are we?
workdir="/data/NC_projects/4120_ARuiz_chicken_variants/analysis"
cd ${workdir}

# create tmp folder
mkdir -p tmpfiles

# how much RAM can java use? - set these to less than your available RAM!
javaopts="-Xms24g -Xmx24g"

#############################################
# GATK4 BAM RECALIBRATION
#############################################

infolder=bam_recalibration
outfolder=HaplotypeCaller
mkdir -p ${outfolder}

dbsnp=reference/ens_dbSNP_150/GRCg6a_dbSNP_150.vcf.gz
reference_fa=reference/Gallus_gallus.GRCg6a.dna.toplevel.fa

# restrict to one chromosome for testing
chr=1

# parallelize the pair hidden Markov models (pair HMM) process
hmmt=12

# run all samples in parallel with 12 threads heach
for bamfile in ${infolder}/5*_mrkdup_srt_recal.bam; do

outpfx=$(basename ${bamfile%_mrkdup_srt_recal.bam})

java ${javaopts} -jar $GATK/gatk.jar \
	HaplotypeCaller  \
	--input ${bamfile} \
	--output ${outfolder}/${outpfx}_${chr}.g.vcf.gz \
	--dbsnp ${dbsnp} \
	--reference ${reference_fa} \
	--intervals ${chr} \
	--emit-ref-confidence GVCF \
	--sample-ploidy 2 \
	--native-pair-hmm-threads ${hmmt} \
	--bam-output ${outfolder}/${outpfx}_${chr}_HC_aligned_reads.bam \
	--tmp-dir tmpfiles/ &

done
