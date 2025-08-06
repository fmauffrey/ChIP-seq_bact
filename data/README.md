# Data Documentation for Snakemake Pipeline

This directory contains the data needed for the analysis. Keep a consistent name for the pipeline to find the correct files.

## Datasets

1. **Reference genome**
   - **Description**: the genome of the reference organism.
   - **Format**: Fasta (.fa)

2. **Annotation**
   - **Description**: annotation of the reference genome.
   - **Format**: GTF (.gtf)

3. **Annotation**
   - **Description**: annotation of the reference genome.
   - **Format**: GFF (.gff)

4. **Experiment data**
   - **Description**: the data files to analyze.
   - **Format**: Fastq (.fq)
   - **Pairing**: add the -1 and -2 tag after the name to refer to forward and reverse, respectively.

## Example

- PA01_GCF_000006765.fa
- PA01_GCF_000006765.gtf
- PA01_GCF_000006765.gff
- pseudo_TF1-1.fq
- pseudo_TF1-2.fq