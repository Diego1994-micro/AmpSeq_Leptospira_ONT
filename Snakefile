###############################################################################
# AMPSEQ ANALYSIS PIPELINE FOR LEPTOSPIRA
# NOTE: This Snakefile is a consensus-calling pipeline for ONT (Nanopore) amplicons, designed to:
# 1. Align reads to a reference genome
# 2. Identify variants (SNPs)
# 3. Determine the consensus base at each position
# 4. Generate a final table of calls by sample and position (REF / ALT / N).
### Interpretation:
### A / C / G / T  → base supported by sufficient evidence
### N              → insufficient coverage or confidence
#
# Useful for comparing samples and determining whether they match or differ from the reference genome.
############################################################################### 


############################################################################### 
# HOW THE SNAKEFILE WORKS
### Each rule defines a step, including its input files, output files, and the command that produces them.
###############################################################################


###############################################################################
# SAMPLE DEFINITION
### Samples must first be defined. A consistent filename prefix (e.g., "barcode" or "sample") facilitates automated batch processing.
# The expression SAMPLES, = glob_wildcards("reads/{sample}.fastq") detects FASTQ filenames and extracts their sample identifiers. A suggested naming convention is sample_Code_Location.fastq.
###############################################################################
SAMPLES, = glob_wildcards("reads/{sample}.fastq")


###############################################################################
# REFERENCE GENOME DEFINITION
### Tara uses the Leptospira interrogans reference genome GCF_000007685.1.
# Another reference genome may be used, provided that the downstream BED file contains coordinates of the regions of interest relative to that genome.
###############################################################################
REF = "reference/GCF_000007685.1_ASM768v1_genomic.fasta"


###############################################################################
# DEFINITION OF INFORMATIVE REGIONS
### The BED file specifies the start and end coordinates of regions containing informative SNPs relative to the L. interrogans GCF_000007685.1 reference genome used by Tara.
# Tara provided this file with all regions analyzed by NAU.
# A different BED file may be used if the regions of interest are mapped to the selected reference genome.
###############################################################################
BED_CALLABLE = 'lepto.new.bed'


###############################################################################
# CONSENSUS BASE-CALLING PARAMETERS
### These variables control how stringent base calling is:
# - A position failing the criteria is assigned "N".
# - Adjusting the thresholds makes consensus calling more or less stringent.

MIN_DP_CALL = 15      # minimum reads required to call a base; otherwise N
MIN_DP_VAR  = 5       # minimum reads to consider a SNP
MIN_ALT_AF  = 0.80    # minimum fraction of ALT-supporting reads (80%)
MIN_ALT_READS = 10    # minimum total ALT-supporting reads
MIN_STRAND_ALT = 3    # minimum ALT-supporting reads per strand
MIN_QUAL = 15         # minimum SNP quality
MIN_MQ = 15           # minimum mapping quality

# ALT -> alternative allele
# REF -> reference allele
# N -> no consensus due to missing information or insufficient sequencing depth
###############################################################################


###############################################################################
# FINAL RULE (rule all)
### The "rule all" target tells Snakemake that analysis is complete when ALL listed files exist.
# In this pipeline, the required outputs are:
### - Per-sample compressed VCF files containing SNPs that passed the specified filters.
### - Corresponding .tbi indexes for those VCF files.
### - Per-sample TSV files containing genomic position, depth, and consensus call, used to build the SNP matrix.
### - The combined consensus-call table and per-sample coverage files.
###############################################################################
rule all:
    input:
        expand("variants/{sample}.consensus.snps.vcf.gz", sample=SAMPLES),
        expand("variants/{sample}.consensus.snps.vcf.gz.tbi", sample=SAMPLES),
        expand("calls/{sample}.consensus_calls.tsv", sample=SAMPLES),
        expand("calls/all_samples_consensus_calls.tsv", sample=SAMPLES),
        expand("coverage/{sample}_coverage.tsv", sample=SAMPLES),
        #"pileup_summary/all_samples_calls.tsv"


###############################################################################
# RULE minimap2_align_sort_index_ont
#
### INPUT
# - FASTQ containing Nanopore reads
# - Reference genome
#
### PROCESS
# Aligns ONT reads to the reference genome and generates a sorted and indexed BAM for each sample.
#
### OUTPUT
# - Sorted and indexed BAM for each sample
###############################################################################
rule minimap2_align_sort_index_ont:
    input:
        reads="reads/{sample}.fastq",
        ref=REF
    output:
        bam="sorted_alignments/{sample}.sorted.bam",
        bai="sorted_alignments/{sample}.sorted.bam.bai"
    threads: 8
    shell:
        r"""
        set -euo pipefail
        mkdir -p sorted_alignments

        minimap2 -ax map-ont -t {threads} {input.ref} {input.reads} \
          | samtools sort -@ {threads} -o {output.bam} -

        samtools index {output.bam}
        """


###############################################################################
# RULE restrict_to_callable_regions
#
### INPUT
# - BED_CALLABLE file specifying AmpSeq regions containing informative SNPs, excluding primer regions
# - Sorted and indexed BAM for each sample
#
### PROCESS
# Retains alignments overlapping the regions specified in the BED file.
#
### OUTPUT
# - Indexed BAM containing alignments overlapping the regions of interest
###############################################################################
rule restrict_to_callable_regions:
    input:
        bam="sorted_alignments/{sample}.sorted.bam",
        bai="sorted_alignments/{sample}.sorted.bam.bai",
        bed=BED_CALLABLE
    output:
        bam="filtered_alignments/{sample}.callable.bam",
        bai="filtered_alignments/{sample}.callable.bam.bai"
    threads: 4
    shell:
        r"""
        set -euo pipefail
        mkdir -p filtered_alignments

        samtools view -@ {threads} -b -L {input.bed} {input.bam} > {output.bam}
        samtools index {output.bam}
        """


###############################################################################
# RULE coverage_per_region
#
### INPUT
# - BED file
# - BAM containing reads in the regions of interest
#
### PROCESS
# Calculates mean coverage per amplicon to identify regions with low coverage.
#
### OUTPUT
# - TSV reporting mean coverage for each AmpSeq target region
###############################################################################
rule coverage_per_region:
    input:
        bam="filtered_alignments/{sample}.callable.bam",
        bed=BED_CALLABLE
    output:
        tsv="coverage/{sample}_coverage.tsv"
    conda: "nanopore"
    shell:
        r"""
        set -euo pipefail
        mkdir -p coverage
        bedtools coverage -a {input.bed} -b {input.bam} -mean > {output.tsv}
        """


###############################################################################
# RULE call_snps_raw
### INPUT
# - BAM with reads in regions of interest
# - Reference genome
#
### PROCESS
# Uses bcftools to identify candidate SNPs without applying the downstream consensus filters.
#
### OUTPUT
# - .vcf.gz: compressed variant-call file containing unfiltered variant candidates
# - .vcf.gz.tbi: Tabix index file
###############################################################################
rule call_snps_raw:
    input:
        bam="filtered_alignments/{sample}.callable.bam",
        bai="filtered_alignments/{sample}.callable.bam.bai",
        ref=REF
    output:
        vcf="variants/{sample}.raw.vcf.gz",
        tbi="variants/{sample}.raw.vcf.gz.tbi"
    conda: "nanopore"
    shell:
        r"""
        set -euo pipefail
        mkdir -p variants

        bcftools mpileup -Ou -f {input.ref} \
          -q 20 -Q 10 \
          -a FORMAT/DP,FORMAT/AD,FORMAT/DP4 \
          {input.bam} \
        | bcftools call --ploidy 1 -m -v -Oz -o {output.vcf}

        bcftools index -t {output.vcf}
        """


###############################################################################
# RULE filter_snps_consensus
### INPUT
# - .vcf.gz: compressed file containing unfiltered variant candidates
# - .vcf.gz.tbi: Tabix index file
#
### PROCESS
# Removes candidate SNPs with inadequate coverage, quality, or alternative-allele fraction.
## The bcftools expression below requires QUAL >= 20, DP >= 5, and ALT allele fraction >= 0.80.
#
### OUTPUT
# - Compressed VCF containing high-confidence SNPs
# - Corresponding Tabix index (.tbi)
###############################################################################
rule filter_snps_consensus:
    input:
        vcf="variants/{sample}.raw.vcf.gz",
        tbi="variants/{sample}.raw.vcf.gz.tbi"
    output:
        vcf="variants/{sample}.consensus.snps.vcf.gz",
        tbi="variants/{sample}.consensus.snps.vcf.gz.tbi"
    conda: "nanopore"
    shell:
        r"""
        set -euo pipefail

        bcftools filter -e '
          TYPE!="snp" ||
          QUAL<20 ||
          FORMAT/DP<5 ||
          (FORMAT/AD[0:1]/FORMAT/DP)<0.8
        ' {input.vcf} -Oz -o {output.vcf}

        bcftools index -t {output.vcf}
        """

###############################################################################
# RULE consensus_call_table
### INPUT
# - VCF containing high-confidence SNPs
# - Tabix index for the VCF
# - Reference genome (REF)
#
### PROCESS
# Builds a table for each callable genomic position, assigning the majority base observed in the reads if it meets the coverage and support thresholds; otherwise assigns N. The VCF is generated independently using the original SNP-calling rules.
#
### OUTPUT
# - One TSV consensus-call table per sample
## Rows: genomic positions with A/C/G/T/N calls
## Columns: Chrom, Pos, Depth, Call
###############################################################################
rule consensus_call_table:
    input:
        bam="filtered_alignments/{sample}.callable.bam",
        bai="filtered_alignments/{sample}.callable.bam.bai",
        ref=REF,
        bed=BED_CALLABLE,
        snps="variants/{sample}.consensus.snps.vcf.gz",
        snps_tbi="variants/{sample}.consensus.snps.vcf.gz.tbi"
    output:
        tsv="calls/{sample}.consensus_calls.tsv"
    run:
        import os
        import subprocess
        import pandas as pd

        os.makedirs("calls", exist_ok=True)

        # ONLY FUNCTIONAL CHANGE: call bases supported by the reads using
        # samtools mpileup. Do not automatically replace positions without
        # an accepted SNP with the reference base.
        # -q filters mapping quality; -Q filters base quality;
        # -aa includes positions with zero coverage.
        cmd = ["samtools", "mpileup", "-aa", "-q", str(MIN_MQ), "-Q", "10",
               "-f", input.ref, "-l", input.bed, input.bam]
        pileup = subprocess.check_output(cmd, text=True)

        def read_bases(pileup_bases, reference):
            """Count A/C/G/T in the pileup; ignore indels and SAMtools markers."""
            counts = {b: 0 for b in "ACGT"}
            i = 0
            while i < len(pileup_bases):
                c = pileup_bases[i]
                if c == "^":
                    i += 2  # '^' followed by a mapping-quality character
                elif c == "$":
                    i += 1
                elif c in "+-":
                    i += 1
                    j = i
                    while i < len(pileup_bases) and pileup_bases[i].isdigit():
                        i += 1
                    if j == i:
                        raise ValueError("Malformed indel in mpileup")
                    i += int(pileup_bases[j:i])
                elif c in ".,":
                    if reference in counts:
                        counts[reference] += 1
                    i += 1
                elif c.upper() in counts:
                    counts[c.upper()] += 1
                    i += 1
                else:
                    # '*', '#', '<', '>', and 'N' are not usable bases
                    i += 1
            return counts

        rows = []
        for line in pileup.splitlines():
            fields = line.split("\t")
            if len(fields) < 5:
                continue
            chrom, pos, refbase = fields[0], int(fields[1]), fields[2].upper()
            counts = read_bases(fields[4], refbase)
            depth = sum(counts.values())
            base, support = max(counts.items(), key=lambda kv: kv[1])
            # Call a base only when supported by sufficient reads.
            # Use N for mixtures or insufficient coverage.
            valid = (depth >= MIN_DP_CALL and support >= MIN_ALT_READS and
                     support / depth >= MIN_ALT_AF)
            rows.append((chrom, pos, depth, base if valid else "N"))

        # Preserve the exact output format of the original Snakefile.
        pd.DataFrame(rows, columns=["Chrom", "Pos", "Depth", "Call"]).to_csv(
            output.tsv, sep="\t", index=False
        )


###############################################################################
# RULE combine_calls
## INPUT
# - TSV files from all previously analyzed samples
#
## PROCESS
# Combines consensus calls from all samples into one matrix.
#
## OUTPUT
# - The main analysis output is "calls/all_samples_consensus_calls.tsv".
## Rows: genomic positions with A/C/G/T/N calls
## Columns: samples
###############################################################################
rule combine_calls:
    input:
        expand("calls/{sample}.consensus_calls.tsv", sample=SAMPLES)
    output:
        "calls/all_samples_consensus_calls.tsv"
    run:
        import pandas as pd
        mats = []
        for f, s in zip(input, SAMPLES):
            df = pd.read_csv(f, sep="\t").set_index(["Chrom","Pos"])
            mats.append(df["Call"].rename(s))
        mat = pd.concat(mats, axis=1).sort_index().fillna("N")
        mat.to_csv(output[0], sep="\t", index=True)
