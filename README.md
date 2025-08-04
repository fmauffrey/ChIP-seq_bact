# Snakemake Pipeline

This project is a Snakemake pipeline designed for data processing. It provides a structured workflow to automate the execution of data analysis tasks.

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