# Config file with run parameters
configfile: "config.yaml"

# Run quality control and trimming pipeline
rule qc:
    message: "Starting quality control and trimming"
    input:
        expand("{sample}/{sample}-1_fastp.fastq", sample=config["input_samples"]),
        expand("{sample}/{sample}-1_fastp_fastqc.zip", sample=config["input_samples"]),

# Reads quality control with Fastp
rule fastp:
    message: "Fastp: {wildcards.sample}"
    input: 
        R1 = "data/{sample}-1.fq",
        R2 = "data/{sample}-2.fq"
    output: 
        R1 = "{sample}/{sample}-1_fastp.fastq",
        R2 = "{sample}/{sample}-2_fastp.fastq",
        html = "02-Fastp/{sample}_fastp.html",
        json = "02-Fastp/{sample}_fastp.json"
    log: "{sample}/{sample}_fastp_log.txt"
    container: "docker://staphb/fastp"
    threads: 1
    params:
        qual = config["fastp"]["qualified_quality_phred"],
        unqual = config["fastp"]["unqualified_percent_limit"],
        ave_qual = config["fastp"]["average_qual"],
        length = config["fastp"]["length_limit"]
    shell:
        "fastp -i {input.R1} -I {input.R2} -o {output.R1} -O {output.R2} -q {params.qual} "
        "-u {params.unqual} -e {params.ave_qual} -l {params.length} -h {output.html} -j {output.json} 2> {log}"

# Reads quality check with FastQC
rule fastqc:
    message: "FastQC: {wildcards.sample}"
    input: 
        "{sample}/{sample}-1_fastp.fastq",
        "{sample}/{sample}-2_fastp.fastq"
    output: 
        "{sample}/{sample}-1_fastp_fastqc.zip",
        "{sample}/{sample}-2_fastp_fastqc.zip"
    container: "docker://staphb/fastqc"
    threads: 1
    shell:
        "fastqc {input} -o {wildcards.sample}"

