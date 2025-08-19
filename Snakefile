# Config file with run parameters
configfile: "config.yaml"

ref_genome = config["reference_genome"]

# Run quality control and trimming pipeline
rule qc:
    message: "Starting quality control and trimming"
    input:
        expand("results/{sample}/QC/{sample}_fastp.fastq", sample=config["input_samples"]),
        expand("results/{sample}/QC/{sample}_fastp_fastqc/fastqc_data.txt", sample=config["input_samples"])

# Rule to align reads on genome
rule align:
    message: "Starting reads alignment"
    input:
        expand("results/{sample}/Align/{sample}.bam", sample=config["input_samples"]),
        expand("results/{sample}/Align/{sample}_phantompeak.txt", sample=config["input_samples"]),

# Rule to call peaks
rule call_peaks:
    message: "Starting peak calling"
    input:
        expand("results/{sample}/Peaks/{sample}_summits.bed", sample=config["input_samples"])

# Rule to annotate peaks
rule annotate_peaks:
    message: "Starting peak annotation"
    input:
        expand("results/{sample}/Annotations/{sample}_closest_genes_summary.txt", sample=config["input_samples"])

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
        min_length = config["fastp"]["min_length"]
    shell:
        "fastp -i {input} -o {output.fastq} -q {params.qual} -u {params.unqual} -e {params.ave_qual} "
        "-l {params.min_length} -h {output.html} -j {output.json} 2> {log}"

# Reads quality check with FastQC
rule fastqc:
    message: "FastQC: {wildcards.sample}"
    input: "results/{sample}/QC/{sample}_fastp.fastq"
    output: "results/{sample}/QC/{sample}_fastp_fastqc/fastqc_data.txt"
    container: "docker://staphb/fastqc"
    threads: 1
    shell:
        "fastqc {input} -o results/{wildcards.sample}/QC --extract --delete"

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

# ChIP-seq alignment quality control
rule phantompeakqualtools:
    message: "Phantompeakqualtools: {wildcards.sample}"
    input: "results/{sample}/Align/{sample}.bam"
    output: "results/{sample}/Align/{sample}_phantompeak.txt"
    container: "docker://seqeralabs/phantompeakqualtools"
    threads: 2
    shell:
        "Rscript scripts/run_spp.R -c={input} -savp -out={output}"

# Peak calling with MACS3
rule macs3:
    message: "MACS3: {wildcards.sample}"
    input: "results/{sample}/Align/{sample}.bam"
    output: 
        peaks="results/{sample}/Peaks/{sample}_peaks.narrowPeak",
        bed="results/{sample}/Peaks/{sample}_summits.bed",
        peaks_xls="results/{sample}/Peaks/{sample}_peaks.xls"
    log: "results/{sample}/Peaks/{sample}_macs3.log"
    threads: 2
    params:
        control=config["control"],
        genome_size=config["macs3"]["genome_size"],
        qvalue=config["macs3"]["qvalue"],
        frag_length=config["macs3"]["fragment_length"]
    shell:
        "macs3 callpeak -t {input} -c results/{params.control}/Align/{params.control}.bam -n {wildcards.sample} "
        "--outdir results/{wildcards.sample}/Peaks --nomodel --extsize {params.frag_length} -f BAM -g {params.genome_size} "
        "-B -q {params.qvalue} 2> {log}"

# Peaks calling summary
rule peaks_summary:
    message: "Generating Peaks calling summary"
    input:
        peaks_QC=expand("results/{sample}/Align/{sample}_phantompeak.txt", sample=config["input_samples"]),
        peaks=expand("results/{sample}/Peaks/{sample}_peaks.narrowPeak", sample=config["input_samples"]),
        bowtie2_log=expand("results/{sample}/Align/{sample}_bowtie2.log", sample=config["input_samples"]),
        fastqc_files=expand("results/{sample}/QC/{sample}_fastp_fastqc/fastqc_data.txt", sample=config["input_samples"])
    output: 
        path="results/peaks_summary.txt"
    threads: 1
    params:
        samples=config["input_samples"],
        genome_size=config["macs3"]["genome_size"]
    script:
        "scripts/peaks_summary.py"

# Rule to find closest genes to peaks
rule bedtools_window:
    message: "Bedtools : {wildcards.sample}"
    input: "results/{sample}/Peaks/{sample}_summits.bed",
    output: "results/{sample}/Annotations/{sample}_closest_genes.txt"
    log: "results/{sample}/Annotations/{sample}_closest_genes.log"
    threads: 1
    container: "docker://staphb/bedtools"
    params:
        genome=config["reference_genome"],
        window=config["bedtools"]["window_size"]
    shell:
        "bedtools window -a {input} -b data/{params.genome}.gff -w {params.window} > {output} 2> {log}"

# Rule to generate a summary of the closest genes
rule closest_genes_summary:
    message: "Closest genes summary: {wildcards.sample}"
    input: 
        closest_genes="results/{sample}/Annotations/{sample}_closest_genes.txt",
        peaks_xls="results/{sample}/Peaks/{sample}_peaks.xls"
    output: 
        path="results/{sample}/Annotations/{sample}_closest_genes_summary.txt"
    threads: 1
    script:
        "scripts/closest_genes_summary.py"
