#!/bin/env bash
# author:Stephane Plaisance (VIB-NC), 2022-03-28

# run deepvariant from docker
# https://github.com/google/deepvariant/blob/r0.9/docs/trio-merge-case-study.md

workdir=/data/NC_projects/4120_ARuiz_chicken_variants/deepvariants

cd ${workdir}

mkdir -p input output tmp_dir

# get copy of input files
#cp -r /data/NC_projects/4120_ARuiz_chicken_variants/gatk4_BP_recal/*_rawmappings_recal.bam input/
#cp -r /data/NC_projects/4120_ARuiz_chicken_variants/gatk4_BP_recal/*_rawmappings_recal.bai input/
#cp /data/biodata/references/galGal6a/Gallus_gallus.GRCg6a.dna.toplevel.fa input/ && \
#samtools faidx input/Gallus_gallus.GRCg6a.dna.toplevel.fa

# set tmp folder until the --intermediate_results_dir bug gets fixed
export TMPDIR="/tmp_dir"

# use almost all threads
nthr=84

# how much RAM can java use? - set these to less than your available RAM!
javaopts="-Xms24g -Xmx24g"

# cycle through 7 samples
for bam in input/*_rawmappings_recal.bam; do

pfx=$(basename ${bam%_rawmappings_recal.bam})

BIN_VERSION="1.3.0"

sudo docker run \
  -v "${workdir}/input":"/input" \
  -v "${workdir}/output":"/output" \
  -v "${workdir}/tmp_dir":"/tmp_dir" \
  google/deepvariant:"${BIN_VERSION}" \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=WGS \
  --ref=/input/Gallus_gallus.GRCg6a.dna.toplevel.fa \
  --reads=/input/"${pfx}"_rawmappings_recal.bam \
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

for gvcf in output/*.g.vcf; do
vcf2index ${gvcf}
done

######################
# merge g.VCF to VCF
######################

myenv=deepvariant
source /etc/profile.d/conda.sh
conda activate ${myenv} || \
  ( echo "# the conda environment ${myenv} was not found on this machine" \
    && exit 1 )

glnexus_cli –config WGS \
  output/*.g.vcf | \
  bcftools view - | \
  bgzip -@ 24 -c \
  > output/merged_7_samples.vcf.gz && \
  tabix -p vcf output/merged_7_samples.vcf.gz

conda deactivate


# 524
c=524
s=0
e=4
d=10
p=20

function getdeepvariantcandidates(){
zcat ${outfolder}/deepvariants_combined_7_snpeff.vcf.gz | \
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

getdeepvariantcandidates
