# Snakemake Pipeline

This project is a Snakemake pipeline designed for the analysis of ChIP-seq data. The pipeline aimed at being simple and straightforward.

## Steps and tools included in the pipeline
### Quality control (FastQC + fastp)
Reads are first filtered using fastp according to the parameters in the config file. A quality control step is then performed using FastQC on the filtered files.

### Reads mapping (Bowtie2, Samtools)

### Peak calling (MACS3) --nomodel is used since macs3 is not able to derive models from these data. A fragment size of around 200 pb is obtained after sonication so --extsize


## Other steps performed manually
- Peak annotation (Homer or ChIPpeakAnno)
- Comparison between samples (Deseq2)
- Motif discovery (RSAT software, peak-motifs command) (MEME) (Homer)
- Comparison between peaks and known TFBS
- Web app for exploring results

## Tips
A too high coverage can lead to many peaks called and this can easily occur with bacteria genomes. Try to subsample in order to get about 1 fold genome coverage, or between 1-10 fold.

## Getting Started

1. **Clone the repository**:
   ```
   git clone <repository-url>
   cd snakemake-pipeline
   ```

2. **Set up the conda environment**:
   ```
   conda env create --name chip-seq --file=envs/environment.yaml
   conda activate chip-seq
   ```

3. **Configure the pipeline**:
   Edit the `config.yaml` file to set your parameters.

   q-value (adjusted with BH) is used for macs3 instead if p-value

4. **Run the reads QC steps**:
   Execute the following command to start QC analysis:
   ```
   snakemake qc --cores <number-of-cores> --use-singularity
   ```
   Pause and check the quality of the filtered fastq files. If OK, proceed to next step. 

5. **Run the alignment steps**:
   Execute the following command to start reads alignment:
   ```
   snakemake align --cores <number-of-cores> --use-singularity
   ```
   Pause and check the quality of the peaks quality (*_phantompeak.txt). If OK, proceed to next step. 

6. **Run the peaks calling steps**:
   Execute the following command to start reads alignment:
   ```
   snakemake call_peaks --cores <number-of-cores> --use-singularity
   ```

## Additional Information

For more details on the data used in this pipeline, please refer to the `data/README.md` file.