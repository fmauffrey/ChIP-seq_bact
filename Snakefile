# Config file with run parameters
configfile: "config.yaml"

# Run the pipeline on all samples listed in the config file
rule run:
    message: "Starting the pipeline"
    input:
        expand("{sample}/{sample}-1_fastqc.zip", sample=config["input_samples"])

# Reads quality check with FastQC
rule fastqc:
    message: "FastQC: {wildcards.sample}"
    input: 
        "data/{sample}-1.fq",
        "data/{sample}-2.fq"
    output: 
        "{sample}/{sample}-1_fastqc.zip",
        "{sample}/{sample}-2_fastqc.zip"
    container: "docker://staphb/fastqc"
    threads: 1
    shell:
        "fastqc {input} -o {wildcards.sample}"