# ChIP-seq analysis pipeline

This project is a Snakemake pipeline designed for the analysis of ChIP-seq data of bacterial strains. The pipeline aimed at being simple and straightforward.

For information about files and format, see [Files and formats](/data/README.md).

## Steps and tools included in the pipeline
### Quality control (FastQC + fastp)
Reads are first filtered using **fastp** according to the parameters in the config file. A quality control step is then performed using **FastQC** on the filtered files.

### Reads mapping (Bowtie2, Samtools)
Filtered reads are mapped against the reference genome using **Bowtie2** aligner. Only erads that mapped once are kept and saved in a .bam alignment file using **Samtools**.

### Peaks calling (MACS3)
Peaks calling is performed with **macs3** in order to detect significant peaks in the genome.

### Peaks annotation (bedtools)
This step identifies all genes located within a specified window around each peak. This approach accounts for variability in transcription start site distances, as the nearest gene is not always the one associated with the annotated transcription site.

## Structure of the config file
- **input_samples** : the list of all the samples to analyze. The control sample must be included in this list.
- **control** : the name (without the extension) of the sample used as control for peak calling.
- **reference_genome** : the name (without the extension) of the reference genome.
- **fastp** : the different parameters for the fastp filtering step.
- **macs3** : the different parameters for the macs3 peak calling step. These parameters should be adapted according to your experiment. In this pipeline, the --nomodel option is used to allow macs3 to call peaks, therefore the --extsize parameter is required.
- **bedtools** : the window parameter for the peak annotation step. 

## Getting Started

1. **Set up the pipeline**:  
   Clone the repository
   ```
   git clone https://github.com/fmauffrey/ChIP-seq_pseudo_zinc
   cd ChIP-seq_pseudo_zinc
   ```
   Copy your files in the data folder then edit the `config.yaml` file to set your parameters.
   ```
   nano config.yaml
   ```

2. **Set up the conda environment**:  
   ```
   conda env create --name chip-seq --file=envs/environment.yaml
   conda activate chip-seq
   ```

3. **Reads QC**:  
   Execute the following command to start QC analysis:
   ```
   snakemake qc --cores <number-of-cores> --use-singularity
   ```
   Pause and check the quality of the filtered fastq files. If OK, proceed to next step.  

4. **Reads alignment on genome**:  
   Execute the following command to start reads alignment:
   ```
   snakemake align --cores <number-of-cores> --use-singularity
   ```
   Pause and check the quality of the peaks quality (*_phantompeak.txt). If OK, proceed to next step. 

5. **Peaks calling**:  
   Execute the following command to start peaks calling:
   ```
   snakemake call_peaks --cores <number-of-cores> --use-singularity
   ```

6. **Metrics summary**:  
   Execute the following command to produce a summary table:
   ```
   snakemake summary --cores <number-of-cores> --use-singularity
   ```

7. **Peaks annotation**  
   Execute the following command to start peaks annotation:
   ```
   snakemake annotate_peaks --cores <number-of-cores> --use-singularity
   ```
   Alternatively, using this command first will execute the entire pipeline from QC to summary file generation.

## Additional information
- A too high coverage can lead to too many peaks called and this can easily occur with bacteria genomes. Try to subsample in order to get about 1-10 fold genome coverage. Use a sub-sampler such as [rasusa](https://github.com/mbhall88/rasusa).