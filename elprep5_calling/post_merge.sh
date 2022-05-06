#!/bin/env bash

#############################################
# GATK4 VARIANT HARD FILTERING
#############################################

# from https://www.nature.com/articles/s41422-020-0349-y
# - MQ < 25.0
# - QUAL < 40.0
# - MQ0 ≥ 4 && ((MQ0/(1.0*DP)) > 0.1)
# --cluster-size 3
# --cluster-window-size 10 (flag more than 3 clustered variants within 10bps +++)

# instructions and filters from:
# https://gatkforums.broadinstitute.org/gatk/discussion/23216/how-to-filter-variants-either-with-vqsr-or-by-hard-filtering
# Genepattern: QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0 || QUAL < 30

# This produces a VCF with the same variant records now annotated with filter status. Specifically, if a record passes all the filters, it receives a PASS label in the FILTER column.
# A record that fails a filter #receives the filter name in the FILTER column, e.g. SOR3.
# If a record fails multiple filters, then each failing filter name appears in the FILTER column separated by semi-colons ; e.g. "MQRankSum-12.5;ReadPosRankSum-8".

infolder=gatk_variants
outfolder=gatk_varianthardfiltering
mkdir -p ${outfolder}

javaopts="-Xms24g -Xmx24g"

###########################################
# 1) hard-Filter SNPs on multiple metrics
###########################################

# produces a VCF with records with SNP-type variants only.
java ${javaopts} -jar $GATK/gatk.jar \
	SelectVariants \
	-V ${infolder}/combined_7.vcf.gz \
	-O ${outfolder}/combined_7_snp.vcf.gz \
	--select-type-to-include SNP

threshold=54.69

java ${javaopts} -jar $GATK/gatk.jar \
	VariantFiltration \
	-V ${outfolder}/combined_7_snp.vcf.gz \
	-O ${outfolder}/combined_7_snp_filtered.vcf.gz \
	--filter "QD < 2.0" --filter-name "QD2" \
	--filter "QUAL < 30.0" --filter-name "QUAL30" \
	--filter "SOR > 3.0" --filter-name "SOR3" \
	--filter "FS > 60.0" --filter-name "FS60" \
	--filter "MQ < 40.0" --filter-name "MQ40" \
	--filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
  	--filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
	--filter-expression "ExcessHet > ${threshold}" --filter-name "ExcessHet" \
    --cluster-size 3 \
    --cluster-window-size 10

######################################################
# 2) hard-Filter INDELs and MIXED on multiple metrics
######################################################

# produces a VCF with records with INDEL-type variants only.
java ${javaopts} -jar $GATK/gatk.jar \
	SelectVariants \
	-V ${infolder}/combined_7.vcf.gz \
	-O ${outfolder}/combined_7_indel.vcf.gz \
	--select-type-to-include INDEL \
	--select-type-to-include MIXED

threshold=54.69

java ${javaopts} -jar $GATK/gatk.jar \
	VariantFiltration \
	-V ${outfolder}/combined_7_indel.vcf.gz \
	-O ${outfolder}/combined_7_indel_filtered.vcf.gz \
	--filter "QD < 2.0" --filter-name "QD2" \
	--filter "QUAL < 30.0" --filter-name "QUAL30" \
	--filter "FS > 200.0" --filter-name "FS200" \
	--filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
	--filter-expression "ExcessHet > ${threshold}" --filter-name "ExcessHet" \
    --cluster-size 3 \
    --cluster-window-size 10

###########################################
## 3) merge SNP and Indel filtered calls
###########################################

# references (indexed)
reference_fa=/data/biodata/references/galGal6a/Gallus_gallus.GRCg6a.dna.toplevel.fa
dbsnp=/data/biodata/references/galGal6a/ens_dbSNP_150/GRCg6a_dbSNP_150.vcf.gz

# combine filtered into one file using Picard
java ${javaopts} -jar $GATK/gatk.jar \
	MergeVcfs \
	-I ${outfolder}/combined_7_snp_filtered.vcf.gz \
	-I ${outfolder}/combined_7_indel_filtered.vcf.gz \
	-R ${reference_fa} \
	-O ${outfolder}/combined_7_snp_indel_filtered.vcf.gz

#############################################
# SnpEff ANNOTATION
#############################################

infolder=${outfolder}
outfolder="snpeff"
mkdir -p ${outfolder}

# remove duplicates introduced when using --select-type-to-include INDEl & MIXED
( zgrep ^"#" ${infolder}/combined_7_snp_indel_filtered.vcf.gz;
zgrep -v ^"#" ${infolder}/combined_7_snp_indel_filtered.vcf.gz \
  | sort -k 1V,1 -k 2n,2 -k 3V,3 -k 4V,4 -k 5V,5 | uniq ) \
  | bgzip -c > ${outfolder}/combined_7_merged.vcf.gz && \
    tabix -p vcf ${outfolder}/combined_7_merged.vcf.gz

build="GRCg6a.105"

java ${javaopts} -jar $SNPEFF/snpEff.jar \
	-htmlStats ${outfolder}/gatk_hardfiltering_snpEff_summary.html \
	-nodownload \
	${build} \
	${outfolder}/combined_7_merged.vcf.gz | \
	bgzip -c > ${outfolder}/combined_7_snpeff.vcf.gz && \
	tabix -p vcf ${outfolder}/combined_7_snpeff.vcf.gz

#############################################
# FILTER AND CREATE CANDIDATE LIST
#############################################

# add hom/het/tot calls for each sample
# count variant in sample1 1:4 only, 5-7 are neutral
#java ${javaopts} -jar $SNPEFF/SnpSift.jar caseControl \
#	"++++000" ${outfolder}/combined_7_snpeff.vcf.gz | \
#	bgzip -c > ${outfolder}/combined_7_snpeff_cnt0.vcf.gz && \
#	tabix -p vcf ${outfolder}/combined_7_snpeff_cnt0.vcf.gz

# same but consider embryos as control samples to get allele counts in 'Controls'
java ${javaopts} -jar $SNPEFF/SnpSift.jar caseControl \
        "++++---" ${outfolder}/combined_7_snpeff.vcf.gz | \
        bgzip -c > ${outfolder}/combined_7_snpeff_cnt.vcf.gz && \
        tabix -p vcf ${outfolder}/combined_7_snpeff_cnt.vcf.gz

# filter PRIVATE calls with the newly created INFO field
zcat ${outfolder}/combined_7_snpeff_cnt.vcf.gz | \
java ${javaopts} -jar $SNPEFF/SnpSift.jar filter \
	"(Cases[0] + Cases[1] == 1)  & \
	((exists ID) & ( ID =~ 'rs' ))" \
	> ${outfolder}/private_candidates.vcf

# create excel table
java ${javaopts} -jar $SNPEFF/SnpSift.jar extractFields \
  -s "," \
  -e "." \
  ${outfolder}/private_candidates.vcf \
  "CHROM" "POS" "ID" "REF" "ALT" "FILTER" \
  "DP" "MQ" \
  "GEN[0].GT" "GEN[1].GT" "GEN[2].GT" "GEN[3].GT" "GEN[4].GT" "GEN[5].GT" "GEN[6].GT" \
  "GEN[0].AD[0]" "GEN[0].AD[1]" "GEN[0].AD[2]" \
  "GEN[1].AD[0]" "GEN[1].AD[1]" "GEN[1].AD[2]" \
  "GEN[2].AD[0]" "GEN[2].AD[1]" "GEN[2].AD[2]" \
  "GEN[3].AD[0]" "GEN[3].AD[1]" "GEN[3].AD[2]" \
  "GEN[4].AD[0]" "GEN[4].AD[1]" "GEN[4].AD[2]" \
  "GEN[5].AD[0]" "GEN[5].AD[1]" "GEN[5].AD[2]" \
  "GEN[6].AD[0]" "GEN[6].AD[1]" "GEN[6].AD[2]" \
  "ANN[0].EFFECT" "ANN[0].GENE" "ANN[0].FEATUREID" "ANN[0].HGVS_C" "ANN[0].HGVS_P" \
  > ${outfolder}/private_candidates.tsv

# also filter Controls
zcat ${outfolder}/combined_7_snpeff_cnt.vcf.gz | \
java ${javaopts} -jar $SNPEFF/SnpSift.jar filter \
        "(Cases[0] + Cases[1] == 1)  & \
         (Controls[0] + Controls[1] < 2) & \
         ((exists ID) & ( ID =~ 'rs' ))" \
        > ${outfolder}/private_candidates_emb.vcf

# create excel table
java ${javaopts} -jar $SNPEFF/SnpSift.jar extractFields \
  -s "," \
  -e "." \
  ${outfolder}/private_candidates_emb.vcf \
  "CHROM" "POS" "ID" "REF" "ALT" "FILTER" \
  "DP" "MQ" \
  "GEN[0].GT" "GEN[1].GT" "GEN[2].GT" "GEN[3].GT" "GEN[4].GT" "GEN[5].GT" "GEN[6].GT" \
  "GEN[0].AD[0]" "GEN[0].AD[1]" "GEN[0].AD[2]" \
  "GEN[1].AD[0]" "GEN[1].AD[1]" "GEN[1].AD[2]" \
  "GEN[2].AD[0]" "GEN[2].AD[1]" "GEN[2].AD[2]" \
  "GEN[3].AD[0]" "GEN[3].AD[1]" "GEN[3].AD[2]" \
  "GEN[4].AD[0]" "GEN[4].AD[1]" "GEN[4].AD[2]" \
  "GEN[5].AD[0]" "GEN[5].AD[1]" "GEN[5].AD[2]" \
  "GEN[6].AD[0]" "GEN[6].AD[1]" "GEN[6].AD[2]" \
  "ANN[0].EFFECT" "ANN[0].GENE" "ANN[0].FEATUREID" "ANN[0].HGVS_C" "ANN[0].HGVS_P" \
  > ${outfolder}/private_candidates_emb.tsv

# DISCVRSeq VariantQC
conda activate java8

reference_fa=/data/biodata/references/galGal6a/Gallus_gallus.GRCg6a.dna.toplevel.fa

java ${javaopts} -jar $BIOTOOLS/DISCVRSeq/DISCVRSeq-1.3.10.jar VariantQC \
  -V ${outfolder}/combined_7_snpeff_cnt.vcf.gz
  -R ${reference_fa} \
  -O ${outfolder}/VariantQC_combined_7_snpeff_cnt.vcf.html

conda deactivate