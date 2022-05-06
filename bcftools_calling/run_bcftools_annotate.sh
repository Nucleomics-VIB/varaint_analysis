#!/bin/env bash
# author:Stephane Plaisance (VIB-NC), 2022-03-28

# http://samtools.github.io/bcftools/howtos/variant-calling.html

workdir=$PWD
cd ${workdir}

outfolder="bcftools_results"
mkdir -p ${outfolder}

###########
# variables
###########

reference_fa="/data/biodata/references/galGal6a/Gallus_gallus.GRCg6a.dna.toplevel.fa"
dbsnp="/data/biodata/references/galGal6a/ens_dbSNP_150/GRCg6a_dbSNP_150.vcf.gz"
chrlist="/data/biodata/references/galGal6a/main_chr.txt"
echo $(seq 1 28; seq 30 33) MT W Z | tr " " "\n" > ${chrlist}
snpeff_build="GRCg6a.105"

javaopts="-Xms24g -Xmx24g"
bcft=4

################
# merge chr data
################

vcfout=bcftools_multi.vcf.gz
bcftools concat --threads ${bcft} \
  ${outfolder}/${outfolder}/*.vcf.gz \
  -O z \
  -o ${outfolder}/${vcfout} && \
  tabix -p vcf ${outfolder}/${vcfout}

######################
# annotate with SnpEff
######################

infile=${vcfout}
outfile=${infile%.vcf.gz}_snpeff.vcf.gz

# the two input files are cases
casecontrolstring="++"

# adding '-nodownload' to use locally built reference
java ${javaopts} -jar $SNPEFF/SnpSift.jar annotate \
  -id ${dbsnp} \
  ${outfolder}/${infile} | \
  java ${javaopts} -jar $SNPEFF/SnpSift.jar \
    caseControl \"${casecontrolstring}\" | \
    java ${javaopts} -jar $SNPEFF/snpEff.jar \
      -htmlStats ${outfile%.vcf.gz}_summary.html \
      -nodownload \
      ${snpeff_build} - | \
        bgzip -c > ${outfolder}/${outfile} && \
        tabix -p vcf ${outfolder}/${outfile}
