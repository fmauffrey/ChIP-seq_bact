# Config file with run parameters
configfile: "config.yaml"

ref_genome = config["reference_genome"]

# Run quality control and trimming pipeline
rule qc:
    message: "Starting quality control and trimming"
    input:
        expand("results/{sample}/QC/{sample}_fastp.fastq", sample=config["input_samples"]),
        expand("results/{sample}/QC/{sample}_fastp_fastqc.zip", sample=config["input_samples"])

# Rule to align reads on genome
rule align:
    message: "Starting the analysis pipeline"
    input:
        expand("results/{sample}/Align/{sample}.bam", sample=config["input_samples"])

# Rule to call peaks
rule call_peaks:
    message: "Starting peak calling"
    input:
        expand("results/{sample}/Peaks/{sample}_peaks.narrowPeak", sample=config["input_samples"])

# Reads quality control with Fastp
rule fastp:
    message: "Fastp: {wildcards.sample}"
    input: "data/{sample}.fastq"
    output: 
        fastq = "results/{sample}/QC/{sample}_fastp.fastq",
        html = "results/{sample}/QC/{sample}_fastp.html",
        json = "results/{sample}/QC/{sample}_fastp.json"
    log: "results/{sample}/QC/{sample}_fastp_log.txt"
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
    input: "results/{sample}/QC/{sample}_fastp.fastq"
    output: "results/{sample}/QC/{sample}_fastp_fastqc.zip"
    container: "docker://staphb/fastqc"
    threads: 1
    shell:
        "fastqc {input} -o results/{wildcards.sample}/QC"

# Indexing reference genome with Bowtie2
rule bowtie2_index:
    message: "Indexing genome with bowtie2"
    input: "data/{ref}.fa".format(ref=ref_genome)
    output: "data/{ref}/{ref}.1.bt2".format(ref=ref_genome)
    container: "docker://staphb/bowtie2"
    threads: 1
    params:
        ref_genome=ref_genome
    shell:
        "bowtie2-build {input} data/{ref_genome}/{ref_genome}"


# Reads mapping with bowtie2
rule bowtie2_align:
    message: "Bowtie2: {wildcards.sample}"
    input: 
        fastq="results/{sample}/QC/{sample}_fastp.fastq",
        index="data/{ref}/{ref}.1.bt2".format(ref=ref_genome)
    output: temp("results/{sample}/Align/{sample}.sam")
    log: "results/{sample}/Align/{sample}_bowtie2.log"
    container: "docker://staphb/bowtie2"
    threads: 2
    params:
        ref_genome=ref_genome
    shell:
        "bowtie2 -x data/{ref_genome}/{ref_genome} -U {input.fastq} "
        "-S {output} --threads {threads} --no-unal --no-mixed --no-discordant 2> {log}"

# Convert SAM to BAM
rule samtools_view:
    message: "Samtools view: {wildcards.sample}"
    input: "results/{sample}/Align/{sample}.sam"
    output: "results/{sample}/Align/{sample}.bam"
    container: "docker://staphb/samtools"
    threads: 2
    shell:
        "samtools view -bS {input} -o {output} && "
        "rm {input}"

# Peak calling with MACS3
rule macs3:
    message: "MACS3: {wildcards.sample}"
    input: "results/{sample}/Align/{sample}.bam"
    output: "results/{sample}/Peaks/{sample}_peaks.narrowPeak"
    log: "results/{sample}/Peaks/{sample}_macs3.log"
    threads: 2
    params:
        control=config["control"],
        genome_size=config["genome_size"],
        qvalue=config["macs3"]["qvalue"]
    shell:
        "macs3 callpeak -t {input} -c results/{params.control}/Align/{params.control}.bam -n {wildcards.sample} "
        "--outdir results/{wildcards.sample}/Peaks --nomodel -f BAM -g {params.genome_size} -B -q {params.qvalue} 2> {log}"