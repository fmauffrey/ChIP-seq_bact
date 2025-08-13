# Config file with run parameters
configfile: "config.yaml"

all_samples = config["input_samples"]
all_samples.append(config["control"])

# Run quality control and trimming pipeline
rule qc:
    message: "Starting quality control and trimming"
    input:
        expand("results/{sample}/{sample}_fastp.fastq", sample=all_samples),
        expand("results/{sample}/{sample}_fastp_fastqc.zip", sample=all_samples)

# Rule to analysis samples
rule analysis:
    message: "Starting the analysis pipeline"
    input:
        expand("results/{sample}/{sample}.sam", sample=all_samples),
        expand("results/{sample}/{sample}_peaks.narrowPeak", sample=config["input_samples"])

# Reads quality control with Fastp
rule fastp:
    message: "Fastp: {wildcards.sample}"
    input: "data/{sample}.fastq"
    output: 
        fastq = "results/{sample}/{sample}_fastp.fastq",
        html = "results/{sample}/{sample}_fastp.html",
        json = "results/{sample}/{sample}_fastp.json"
    log: "results/{sample}/{sample}_fastp_log.txt"
    container: "docker://staphb/fastp"
    threads: 1
    params:
        qual = config["fastp"]["qualified_quality_phred"],
        unqual = config["fastp"]["unqualified_percent_limit"],
        ave_qual = config["fastp"]["average_qual"],
        length = config["fastp"]["length_limit"]
    shell:
        "fastp -i {input} -o {output.fastq} -q {params.qual} -u {params.unqual} -e {params.ave_qual} "
        "-l {params.length} -h {output.html} -j {output.json} 2> {log}"

# Reads quality check with FastQC
rule fastqc:
    message: "FastQC: {wildcards.sample}"
    input: "results/{sample}/{sample}_fastp.fastq"
    output: "results/{sample}/{sample}_fastp_fastqc.zip"
    container: "docker://staphb/fastqc"
    threads: 1
    shell:
        "fastqc {input} -o results/{wildcards.sample}"

# Indexing reference genome with Bowtie2
rule bowtie2_index:
    message: "Indexing genome with bowtie2"
    input: "data/{params.ref_genome}.fa"
    output: "data/{params.ref_genome}/{params.ref_genome}.1.bt2"
    container: "docker://staphb/bowtie2"
    threads: 1
    params:
        ref_genome = config["reference_genome"]
    shell:
        "bowtie2-build {input} data/{params.ref_genome}/{params.ref_genome}"


# Reads mapping with bowtie2
rule bowtie2_align:
    message: "Bowtie2: {wildcards.sample}"
    input: 
        fastq="results/{sample}/{sample}_fastp.fastq",
        index="data/{params.ref_genome}/{params.ref_genome}.1.bt2"
    output: "results/{sample}/{sample}.sam"
    log: "results/{sample}/{sample}_bowtie2.log"
    container: "docker://staphb/bowtie2"
    threads: 2
    params:
        ref_genome = config["reference_genome"]
    shell:
        "bowtie2 -x data/{params.ref_genome}/{params.ref_genome} -1 {input.fastq} "
        "-S {output} --threads {threads} --no-unal --no-mixed --no-discordant 2> {log}"

# Peak calling with MACS3
rule macs3:
    message: "MACS3: {wildcards.sample}"
    input: "results/{sample}/{sample}.sam"
    output: "results/{sample}/{sample}_peaks.narrowPeak"
    log: "results/{sample}/{sample}_macs3.log"
    threads: 2
    shell:
        "macs3 callpeak -t {input} -f SAM -g hs -n {wildcards.sample} --outdir {wildcards.sample} "
        "--nomodel --shift -100 --extsize 200 --keep-dup all 2> {log}"