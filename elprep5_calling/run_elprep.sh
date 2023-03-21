#!/bin/bash

# script: run_elprep.sh

# conda activate elprep
# adapt next line to point to the right conda.sh init script
# see conda activate script for details
source /etc/profile.d/conda.sh
conda activate elprep || \
  ( echo "# the conda environment 'elprep' was not found on this machine" ;
    echo "# please read the top part of the script!" \
    && exit 1 )

# IO
reffolder=/data/biodata/references/galGal6a
ref=${reffolder}/Gallus_gallus.GRCg6a.dna.toplevel.elfasta
dbsnp=${reffolder}/ens_dbSNP_150/GRCg6a_dbSNP_150.elsites

pixeldist=2500

# threads, adapt to your machine
nthr=22

workdir=/data/NC_projects/4120_ARuiz_chicken_variants
cd ${workdir}
mkdir -p tmpfiles elprep

######################
# Memory and TMP Files
######################

# handle tmpfiles as output for intermediate files if sfm mode
#elprepmode=filter # all in RAM
elprepmode=sfm # less in RAM

if [ ${elprepmode} == "sfm" ]; then
tmpopt="--tmp-path ${workdir}/tmpfiles"
else
tmpopt=""
fi

# run loop
for inbam in bwa_mappings/*_rawmappings.bam; do

pfx=$(basename ${inbam%_rawmappings.bam})
outbam=${workdir}/elprep/${pfx}_elprep_mappings.bam

elprep ${elprepmode} ${inbam} ${outbam} \
  --optical-duplicates-pixel-distance ${pixeldist} \
  --mark-duplicates \
  --mark-optical-duplicates ${workdir}/elprep/${pfx}.output.metrics \
  --sorting-order coordinate \
  --bqsr ${workdir}/elprep/${pfx}_bqsr.output.recal \
  --known-sites ${dbsnp} \
  --reference ${ref}  \
  --haplotypecaller ${workdir}/elprep/${pfx}.g.vcf.gz \
  --reference-confidence GVCF \
  --nr-of-threads ${nthr} \
  --timed \
  --log-path ${workdir}/elprep \
  ${tmpopt}

done
