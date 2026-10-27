AmpSeq Nanopore Consensus Pipeline

This Snakemake pipeline analyzes Oxford Nanopore (ONT) AmpSeq reads from Leptospira samples. It aligns reads to a reference genome, identifies SNPs, and produces consensus base calls for selected genomic regions.

Requirements
- Snakemake, minimap2, SAMtools, BCFtools, BEDTools, and pandas
- A working nanopore environment for rules that specify conda: "nanopore"
- FASTQ reads in reads/ (one file per sample, e.g., sample_01.fastq)
- Reference FASTA at reference/GCF_000007685.1_ASM768v1_genomic.fasta
- BED file lepto.new.bed with target-region coordinates matching that reference genome (primer regions excluded)

Adjust the input paths in the Snakefile if your files have different names or locations.


Run

From the project directory, with the required software and environment available:
snakemake -s Snakefile --cores 8
Rename the provided English Snakefile to Snakefile, or replace Snakefile in the command with its actual filename.


Workflow and outputs
1. Alignment: minimap2 maps reads to the reference; samtools sorts and indexes BAM files (sorted_alignments/).
2. Region selection and coverage: Retain alignments overlapping BED targets (filtered_alignments/) and calculate mean coverage per region (coverage/).
3. SNP calling: bcftools generates raw VCFs and filters SNPs by quality, depth, and alternative-allele fraction (variants/).
4. Consensus calls: samtools mpileup counts observed bases at each targeted position. A base is assigned only with at least 15 usable reads in total, 10 reads supporting the most frequent base, and 80% agreement; otherwise, the call is N (calls/).
5. Combined matrix: calls/all_samples_consensus_calls.tsv combines the consensus calls across samples for downstream comparisons and phylogenetic analyses.


Important: Positions without an accepted SNP are not automatically filled with the reference base. Consensus bases are determined from the reads themselves. SNP filtering and consensus calling are separate steps with their own thresholds.
Reference coordinates

The supplied configuration uses L. interrogans GCF_000007685.1. If another reference genome is used, the BED coordinates must be updated accordingly.
