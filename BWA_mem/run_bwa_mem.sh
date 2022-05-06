#!/bin/env bash

# map paired-reads to reference
# author:Stephane Plaisance (VIB-NC), 2022-03-18
# from: https://github.com/BITS-VIB/NGS-Variant-Analysis-training-2020

workdir=$PWD
cd ${workdir}

# multithreading, adapt to your cpu count
bwathr=84
samtoolsthr=24

# reference genome (.fai and .dict accessory files are present in the same folder)
reference_fa=/data/biodata/references/galGal6a/Gallus_gallus.GRCg6a.dna.toplevel.fa

# a BWA index exists for the reference genome (else create it with bwa index)
bwaidxfolder=/data/biodata/references/galGal6a
bwaidx=GRCg6a.105

#############################################
# BWA mem MAP read pairs
#############################################

outfolder=${workdir}/bwa_mappings
mkdir -p ${outfolder}

# samblaster required to mark duplicates
myenv=samblaster
source /etc/profile.d/conda.sh
conda activate ${myenv} || \
  ( echo "# the conda environment ${myenv} was not found on this machine" \
    && exit 1 )

# loop through read pairs
for reads_1 in reads/*_R1_001.fastq.gz; do

reads_2=${reads_1%_R1_001.fastq.gz}_R2_001.fastq.gz
outpfx=$(basename ${reads_1%_R1_001.fastq.gz})

# map using BWA mem and mark PCR duplicate reads (run only once)
if [ ! -f ${outfolder}/${outpfx}_mapping_done ]; then
	echo "# mapping ${outpfx} reads with BWA mem"
	RSTRING="'@RG\tID:"${outpfx}"\tSM:"${outpfx}"'"
	
	cmd="bwa mem -t ${bwathr} -M \
	  -R ${RSTRING} \
	  ${bwaidxfolder}/${bwaidx} \
	  ${reads_1} \
	  ${reads_2} \ 
	    | samblaster -M \
	    | samtools view -Sb - -o ${outfolder}/${outpfx}_mkrdup.bam"
	echo "# ${cmd}"
	eval ${cmd} 
	
	# sort and index
	samtools sort -@ ${samtoolsthr} \
	  ${outfolder}/${outpfx}_mkrdup.bam \
	  -o ${outfolder}/${outpfx}_mappings.bam && \
	    samtoolls index ${outfolder}/${outpfx}_mappings.bam

	# get flagstats
	# then write 'mapping_done' flag
	samtools flagstat -@ ${samtoolsthr} \
		${outfolder}/${outpfx}_mappings.bam \
		> ${outfolder}/${outpfx}_mappings_flagstats.txt && \
		touch ${outfolder}/${outpfx}_mapping_done
else
	echo "# BWA mapping already done for ${outpfx}"
fi

done
