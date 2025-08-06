# Snakemake Pipeline

This project is a Snakemake pipeline designed for the analysis of ChIP-seq data. The pipeline aimed at being simple and straightforward.

## Steps and tools
- Quality control (FastQC + fastp)
- Reads mapping (Bowtie2)
- Peak calling (MACS3) Careful parameters optimization is required
- Peak annotation (Homer)
- Comparison between samples (Deseq2)

## Other possible steps


## Tips
A too high coverage can lead to many peaks called and this can easily occur with bacteria genomes. Try to subsample in order to get about 1 fold genome coverage, or between 1-10 fold.

## Project Structure

- **Snakefile**: Defines the workflow, including rules, input/output files, and commands for each step.
- **config.yaml**: Contains configuration parameters for the pipeline, allowing users to define variables used throughout the workflow.
- **envs/environment.yaml**: Specifies the conda environment with required packages and their versions for consistent execution.
- **scripts/process_data.py**: Implements the data processing logic, called by Snakemake rules to perform specific tasks.
- **data/README.md**: Documentation about the datasets used in the pipeline, including descriptions and formats.

## Getting Started

1. **Clone the repository**:
   ```
   git clone <repository-url>
   cd snakemake-pipeline
   ```

2. **Set up the conda environment**:
   ```
   conda env create -f envs/environment.yaml
   conda activate <environment-name>
   ```

3. **Configure the pipeline**:
   Edit the `config.yaml` file to set your parameters.

4. **Run the pipeline**:
   Execute the following command to start the Snakemake workflow:
   ```
   snakemake --cores <number-of-cores>
   ```

## Additional Information

For more details on the data used in this pipeline, please refer to the `data/README.md` file.