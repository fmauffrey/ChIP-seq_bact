# Files and formats

This directory contains the data needed for the analysis. Put all files (samples, control, genome) in this folder. Keep a consistent name for the pipeline to find the correct files. A control sample is expected for noise reduction.

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

4. **Sample data**
   - **Description**: the data files to analyze and the control sample. The pipeline has been designed for single-End sequencing data.
   - **Format**: Fastq (.fastq)

## Example

- PA01_GCF_000006765.fa
- PA01_GCF_000006765.gtf
- PA01_GCF_000006765.gff
- pseudo_TF1.fastq
- control_TF1.fastq