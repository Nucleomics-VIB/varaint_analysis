#!/bin/env bash
# author:Stephane Plaisance (VIB-NC), 2019-12-12
# from: https://github.com/BITS-VIB/NGS-Variant-Analysis-training-2020

## Where are we?
workdir="/data/NC_projects/4120_ARuiz_chicken_variants/analysis"
cd ${workdir}

# create tmp folder
mkdir -p tmpfiles
export TMPDIR=tmpfiles

# multithreading, adapt to your cpu count
bwathr=84
samtoolsthr=24

# how much RAM can java use? - set these to less than your available RAM!
javaopts="-Xms24g -Xmx24g"

# samblaster required to mark duplicates
myenv=samblaster
source /etc/profile.d/conda.sh
conda activate ${myenv} || \
  ( echo "# the conda environment ${myenv} was not found on this machine" \
    && exit 1 )

#############################################
# BWA index & MEM mapping
#############################################

## create a bwa index when absent
bwaidxfolder=reference
bwaidx="${bwaidxfolder}/GRCg6a.105"

## map reads to reference

for reads_1 in reads/*_R1.fastq.gz; do

reads_2=${reads_1%_R1.fastq.gz}_R2.fastq.gz

# edit in the command below
outpfx=$(basename ${reads_1%_R1.fastq.gz})

outfolder=mappings
mkdir -p ${workdir}/${outfolder}

RG=\''@RG\tID:'${outpfx}'\tSM:'${outpfx}\'

# version using samblaster for marking duplicates
cmd="bwa mem \
	-t ${bwathr} \
	-M \
	-R ${RG} \
	${bwaidx} \
	${reads_1} \
	${reads_2} | \
	samblaster -M | \
	${samtools} view -Sb - -o ${outfolder}/${outpfx}_dupmrk.bam"

echo "# ${cmd}"
eval ${cmd}

samtools sort -@ ${samtoolsthr} ${outfolder}/${outpfx}_dupmrk.bam \
  > ${outfolder}/${outpfx}_mappings_dupmrk.bam && \
  samtools index ${outfolder}/${outpfx}_mappings_dupmrk.bam

samtools flagstat -@ ${samtoolsthr} \
		${outfolder}/${outpfx}_mappings_dupmrk.bam \
		> ${outfolder}/${outpfx}_mappings_dupmrk_flagstats.txt

done

# cleanup
conda deactivate
