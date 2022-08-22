# variant_analysis
basic scripts to perform small variant analysis on a server

# Read mapping

## BWA mem for read mapping

# Variant calling

## Several alternatives for variant calling

We aimed at producing gVCF data or call all genomes simultaneously in order to handle Ref-calls (same as reference genome) and No-calls (no read support) and obtain a better multigenome dataset. It is also possible to call each genome to VCF and merge the VCF files but doing so the No-calls will be represented as Ref-calls **which is not good** and is what we get below with ```freebayes-parallel```.

### classical bcftools calling

NOTE: We use here ```bcftools mpileup``` (with new arguments) rather than the old ```samtools mpileup``` syntax used in older tutorials.

### GATK4 best ptractice (very slow single core steps)

### elPrep5 GATK4 accellerated alternative (developped for human genomes)

### Google deepvariant alternative (very popular)

### freebayes alternatrive (parallel version used here for speedup)

NOTE: Freebayes is used here as it is a very poular tool. However freebayes does not (yet) seem to produce gVCF to a level comparable to other tools, we are therefore producing separate VCF files here and merging them with bcftools. A particular issue is when a variant in one genome is a SNV while an indel is present in another genome. In such case, merging will not reflect the complexity of the situation and might lead to call issues.

# Variant annotation and filtering

## SNPEFF for variant annotation and SNPSIFT for filtering

