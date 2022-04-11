# variant_analysis
basic scripts to perform small variant analysis on a server

## BWA mem for read mapping

## Several alternatives for variant calling

We aimed at producing gVCF data or call all genomes simultaneously in order to handle Ref-calls (same as reference genome) and No-calls (no read support) and obtain a better multigenome dataset. It is also possible to call each genome to VCF and merge the VCF files but doing so the No-calls will be represented as Ref-calls which is not good.

### classical bcftools calling

NOTE: We use here ```bcftools mpileup``` (with new arguments) rather than the old ```samtools mpileup``` syntax used in older tutorials.

### GATK4 best ptractice (very slow single core steps)

### elPrep5 GATK4 alternative (developped for human genomes)

### freebayes alternatrive (parallel version used heer for speedup)

### Google deepvariant alternative (very popular)

## SNPEFF for variant annotation and SNPSIFT for filtering

