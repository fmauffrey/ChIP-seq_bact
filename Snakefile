# Config file with run parameters
configfile: "config.yaml"

# Run quality control and trimming pipeline
rule qc:
    message: "Starting quality control and trimming"
    input:
        expand("{sample}/{sample}-1_fastp.fastq", sample=config["input_samples"]+config["control"]),
        expand("{sample}/{sample}-1_fastp_fastqc.zip", sample=config["input_samples"]+config["control"])

# Rule to analysis samples
rule analysis:
    message: "Starting the analysis pipeline"
    input:
        expand("{sample}/{sample}.sam", sample=config["input_samples"]+config["control"]),
        expand("{sample}/{sample}_peaks.narrowPeak", sample=config["input_samples"])

# Reads quality control with Fastp
rule fastp:
    message: "Fastp: {wildcards.sample}"
    input: 
        R1 = "data/{sample}-1.fq",
        R2 = "data/{sample}-2.fq"
    output: 
        R1 = "{sample}/{sample}-1_fastp.fastq",
        R2 = "{sample}/{sample}-2_fastp.fastq",
        html = "{sample}/{sample}_fastp.html",
        json = "{sample}/{sample}_fastp.json"
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

# Indexing reference genome with Bowtie2
rule bowtie2_index:
    message: "Indexing genome with bowtie2"
    input: f"data/{config['reference_genome']}.fa"
    output: f"data/{config['reference_genome']}/{config['reference_genome']}.1.bt2"
    container: "docker://staphb/bowtie2"
    threads: 1
    shell:
        f"bowtie2-build {input} data/{config['reference_genome']}/{config['reference_genome']}"


# Reads mapping with bowtie2
rule bowtie2_align:
    message: "Bowtie2: {wildcards.sample}"
    input: 
        R1="{sample}/{sample}-1_fastp.fastq",
        R2="{sample}/{sample}-2_fastp.fastq",
        index=f"data/{config['reference_genome']}/{config['reference_genome']}.1.bt2"
    output: 
        "{sample}/{sample}.sam"
    log: "{sample}/{sample}_bowtie2.log"
    container: "docker://staphb/bowtie2"
    threads: 2
    shell:
        "bowtie2 -x data/{config[reference_genome]}/{config[reference_genome]} -1 {input.R1} -2 {input.R2} "
        "-S {output} --threads {threads} --no-unal --no-mixed --no-discordant "
        "2> {log}"

# Peak calling with MACS3
rule macs3:
    message: "MACS3: {wildcards.sample}"
    input: 
        "{sample}/{sample}.sam"
    output: 
        "{sample}/{sample}_peaks.narrowPeak"
    log: "{sample}/{sample}_macs3.log"
    threads: 2
    shell:
        "macs3 callpeak -t {input} -f SAM -g hs -n {wildcards.sample} --outdir {wildcards.sample} "
        "--nomodel --shift -100 --extsize 200 --keep-dup all 2> {log}"