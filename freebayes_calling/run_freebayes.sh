#!/bin/env bash
# author:Stephane Plaisance (VIB-NC), 2022-03-28

# run freebayes from conda (parallel)

workdir="/data/NC_projects/4120_ARuiz_chicken_variants"

javaopts="-Xms24g -Xmx24g"
nthr=84
bcft=4

reference_fa="/data/biodata/references/galGal6a/Gallus_gallus.GRCg6a.dna.toplevel.fa"
dbsnp="/data/biodata/references/galGal6a/ens_dbSNP_150/GRCg6a_dbSNP_150.vcf.gz"

cd ${workdir}

# create tmp folder
mkdir -p tmpfiles

####################
## Freebayes caller
####################

infolder=${workdir}/gatk4_BP_recal
outfolder=${workdir}/freebayes_vcf/individual_calls
mkdir -p ${outfolder}

# conda env with freebayes and bamtools installed (+ samtools & bcftools)
myenv=freebayes
source /etc/profile.d/conda.sh
conda activate ${myenv} || \
  ( echo "# the conda environment ${myenv} was not found on this machine" \
    && exit 1 )

for bam in ${infolder/*_rawmappings_recal.bam}; do

pfx=$(basename ${bam%_rawmappings_recal.bam})

# create coverage-weighted regions
bamtools coverage -in ${bam} \
  | coverage_to_regions.py ${reference_fa}.fai 500 \
  > tmpfiles/${reference_fa}.500.regions

# call variants using multiple cores
freebayes-parallel tmpfiles/${reference_fa}.500.regions ${nthr} \
  -f ${reference_fa} ${bam} \
  > ${outfolder}/${pfx}_freebayes.vcf

done

conda deactivate

########################
## Freebayes normalize
########################

# normalize freebayes calls
infolder=${workdir}/freebayes_vcf/individual_calls
outfolder=${workdir}/freebayes_vcf/bcftools_merged_norm-calls
mkdir -p ${outfolder}

for vcf in ${infolder}/*_freebayes.cf.gz; do

pfx=$(basename ${vcf%.vcf.gz})

# normalization also produces Ref only rows
# filter out such rows before merging

bcftools norm \
  --threads ${bcft} \
  -f ${reference_fa} \
  ${vcf} | \
  java ${javaopts} -jar $SNPEFF/SnpSift.jar filter "countVariant()>0" | \
  bgzip -c > ${outfolder}/${pfx}_norm.vcf.gz && \
  tabix -p vcf ${outfolder}/${pfx}_norm.vcf.gz

done

####################
## Freebayes merge
####################

bcftools merge \
  --threads ${bcft} \
  -O v \
  -o >(cat | bgzip -c > ${outfolder}/freebayes_merged.vcf.gz) \
  ${infolder}/524_freebayes_norm.vcf.gz \
  ${infolder}/526_freebayes_norm.vcf.gz \
  ${infolder}/528_freebayes_norm.vcf.gz \
  ${infolder}/530_freebayes_norm.vcf.gz \
  ${infolder}/W201120E14-EMB_freebayes_norm.vcf.gz \
  ${infolder}/W291020E14-EMB_freebayes_norm.vcf.gz \
  ${infolder}/W291020E15-EMB_freebayes_norm.vcf.gz && \
  tabix -p vcf ${outfolder}/freebayes_merged.vcf.gz

####################
## identify dbSNP
####################

# add dbSNP annotations
# java ${javaopts} -jar $SNPEFF/SnpSift.jar annotate \
#   -id ${dbsnp} \
#   ${outfolder}/freebayes_merged.vcf.gz | \
#   bgzip -c > ${outfolder}/freebayes_merged_ID.vcf.gz && \
#   tabix -p vcf ${outfolder}/freebayes_merged_ID.vcf.gz

##########################################################
# ADD GROUP ANNOTATION, FILTER AND CREATE CANDIDATE LIST
##########################################################

# java ${javaopts} -jar $SNPEFF/SnpSift.jar caseControl \
#         "++++---" ${outfolder}/freebayes_merged_ID.vcf.gz | \
#         bgzip -c > ${outfolder}/freebayes_merged_ID_cnt.vcf.gz && \
#         tabix -p vcf ${outfolder}/freebayes_merged_ID_cnt.vcf.gz
  
#########################
# ADD SnpEff ANNOTATION
#########################

# build="GRCg6a.105"
# 
# java ${javaopts} -jar $SNPEFF/snpEff.jar \
# 	-htmlStats ${outfolder}/freebayes_snpEff_summary.html \
# 	-nodownload \
# 	${build} \
# 	${outfolder}/freebayes_merged_ID_cnt.vcf.gz | \
# 	bgzip -c > ${outfolder}/freebayes_combined_7_snpeff.vcf.gz && \
# 	tabix -p vcf ${outfolder}/freebayes_combined_7_snpeff.vcf.gz


##############################################################
## identify dbSNP + ADD GROUP ANNOTATION + SnpEff ANNOTATION
##############################################################

build="GRCg6a.105"
infile=${outfolder}/freebayes_merged.vcf.gz
outfile=${outfolder}/freebayes_combined_7_snpeff.vcf.gz

java ${javaopts} -jar $SNPEFF/SnpSift.jar annotate \
  -id ${dbsnp} \
  ${infile} | \
  java ${javaopts} -jar $SNPEFF/SnpSift.jar caseControl "++++---" | \
  java ${javaopts} -jar $SNPEFF/snpEff.jar \
	-htmlStats ${outfile%.vcf.gz}_summary.html \
	-nodownload \
	${build} - | \
	bgzip -c > ${outfile} && \
	tabix -p vcf ${outfile}


###########################
# FILTER candidates
###########################

# QUAL > 

# dbSNPID exists

# 1 sample is Variant
# sample alt depth > 10
# sample GL <= -3

# 0 or 1 embr is variant
# embr alt depth > 0
# embr GL <= -3

# 524
c=524
s=0
e=4
d=10
p=-3

function getcandidates(){
zcat ${outfolder}/freebayes_combined_7_snpeff.vcf.gz | \
java ${javaopts} -jar $SNPEFF/SnpSift.jar filter \
        "(( na FILTER ) | (FILTER = 'PASS')) & \
         ((exists ID) & ( ID =~ 'rs' )) & \
         (Cases[0] + Cases[1] == 1) & \
	     (Controls[0] + Controls[1] < 2) & \
         (GEN[${s}].AO > ${d}) & \
         (GEN[${s}].GL[0] <= ${p}) & \
         (GEN[${s}].GL[1] <= ${p}) & \
         (GEN[${e}].AO > 0) & \
         (GEN[${e}].GL[0] <= ${p}) & \
         (GEN[${e}].GL[1] <= ${p})" \
         > ${outfolder}/freebayes_filtered_${c}.vcf
}



getcandidates
# => 229
 
# 526
c=526
s=1
e=5
d=10
p=-3

getcandidates
# => 205

# 528
c=528
s=2
e=6
d=10
p=-3

getcandidates
# => 233

# 530 (no paired embryo)
c=530
s=3
d=20
p=-6
zcat ${outfolder}/freebayes_combined_7_snpeff.vcf.gz | \
java ${javaopts} -jar $SNPEFF/SnpSift.jar filter \
        "(( na FILTER ) | (FILTER = 'PASS')) & \
         ((exists ID) & ( ID =~ 'rs' )) & \
         (Cases[0] + Cases[1] == 1) & \
	     (Controls[0] + Controls[1] == 0) & \
         (GEN[${s}].AO > ${d}) & \
         (GEN[${s}].GL[0] <= ${p}) & \
         (GEN[${s}].GL[1] <= ${p})" \
         > ${outfolder}/freebayes_filtered_${c}.vcf

# => 12013


#########################################
# RUN DISCVRSeq VariantQC
#########################################

myenv=java8
source /etc/profile.d/conda.sh
conda activate ${myenv} || \
  ( echo "# the conda environment ${myenv} was not found on this machine" \
    && exit 1 )

reference_fa=/data/biodata/references/galGal6a/Gallus_gallus.GRCg6a.dna.toplevel.fa

java ${javaopts} -jar $BIOTOOLS/DISCVRSeq/DISCVRSeq-1.3.10.jar VariantQC \
  -V ${outfolder}/freebayes_merged_ID_cnt.vcf.gz \
  -R ${reference_fa} \
  -O ${outfolder}/VariantQC_freebayes_combined_7_snpeff.vcf.html

conda deactivate

