#!/bin/env bash

# script: run_deepvariant_1.5.0.sh
# author:Stephane Plaisance (VIB-NC), 2023-03-19

# run deepvariant from docker v1.5.0
# get with: docker pull google/deepvariant
# commands from trio tutorial
# https://github.com/google/deepvariant/blob/r0.9/docs/trio-merge-case-study.md

workdir="/data/NC_projects/4395_ARuiz_chicken_variants"
cd ${workdir}

# set tmp folder until the --intermediate_results_dir bug gets fixed
export tmpdir="/tmp_dir"
mkdir -p ${tmpdir}

# use almost all threads
nthr=84

bamfolder="bwa_mappings"
reffolder="reference"
outfolder="deepvariant_results"
mkdir -p ${outfolder}

# samtools faidx indexed fasta
reference="Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa"

# Replace this string with exactly one of the following [WGS,WES,PACBIO,ONT_R104,HYBRID_PACBIO_ILLUMINA]
model_type="WGS"

# cycle through 7 samples
for bamfile in ${bamfolder}$/*_rawmappings_recal.bam; do

bam=$(basename ${bamfile})
pfx=${bam%_rawmappings_recal.bam}

BIN_VERSION="1.5.0"

sudo docker run \
  -v "${workdir}/${BAMFOLDER}":"/inbam" \
  -v "${workdir}/${REFFOLDER}":"/inref" \
  -v "${workdir}/${OUTFOLDER}":"/output" \
  -v "${workdir}/${TMPDIR}":"/tmp_dir" \
  google/deepvariant:"${BIN_VERSION}" \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=${model_type} \
  --ref=/inref/${reference} \
  --reads=/inbam/${bam} \
  --output_vcf=/output/"${pfx}".vcf \
  --output_gvcf=/output/"${pfx}".g.vcf \
  --num_shards="${nthr}" \
  --intermediate_results_dir /tmp_dir \
  --logging_dir=/output/"${pfx}"_logs \
  --dry_run=false

done

# stop here, the remaining will be run later
exit 0

######################
# compress and index
######################

for gvcf in ${outfolder}/*.g.vcf; do
vcf2index ${gvcf}
done

##########################################
# merge g.VCF to VCF using GLnexus 1.4.1 #
##########################################

myenv=deepvariant
source /etc/profile.d/conda.sh
conda activate ${myenv} || \
  ( echo "# the conda environment ${myenv} was not found on this machine" \
    && exit 1 )

# install libjemalloc-dev (ubuntu)
# preload installed library at runtime

LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libjemalloc.so glnexus_cli \
  –config WGS \
  ${outfolder}/*.g.vcf | \
  bcftools view - | \
  bgzip -@ 24 -c \
  > ${outfolder}/merged_samples.vcf.gz && \
  tabix -p vcf output/merged_samples.vcf.gz

conda deactivate

exit 0

################################
# filter variants with SnpSift #
################################

# 524
c=524
s=0
e=4
d=10
p=20

# how much RAM can java use? - set these to less than your available RAM!
javaopts="-Xms24g -Xmx24g"

function getdeepvariantcandidates(){
zcat ${outfolder}/merged_samples.vcf.gz | \
java ${javaopts} -jar $SNPEFF/SnpSift.jar filter \
        "(( na FILTER ) | (FILTER = 'PASS')) & \
         ((exists ID) & ( ID =~ 'rs' )) & \
         (Cases[0] + Cases[1] == 1) & \
         (Controls[0] + Controls[1] < 2) & \
         (GEN[${s}].AD[1] > ${d}) & \
         (GEN[${s}].PL[0] >= ${p}) & \
         (GEN[${s}].PL[1] >= ${p}) & \
         (GEN[${e}].AD[1] > 0)" \
         > ${outfolder}/deepvariants_filtered_${c}.vcf
}

# run function here
getdeepvariantcandidates
