import pandas as pd
import re


def add_mapping_reads(summary):
    # Extract the number of reads mapped to the genome
    summary["Mapping reads"] = []
    for bowtie2_log_file in snakemake.input.bowtie2_log:
        with open(bowtie2_log_file, "r") as f:
            text = f.read()
            value = re.search(r'(\d+)\s*\([^)]+\)\s+aligned exactly 1 time', text).group() # To improve
            reads_number = re.search(r'^\d+', value).group()
            summary["Mapping reads"].append(reads_number)
    return summary

def add_phantompeak(summary):
    # Extract alignments QC metrics and add them to summary
    summary["NSC"] = []
    summary["RSC"] = []
    summary["QualityTag"] = []
    for QC_file in snakemake.input.peaks_QC:
        with open(QC_file, "r") as f:
            line = f.read().strip("\n").split("\t")
            summary["NSC"].append(line[8])
            summary["RSC"].append(line[9])
            summary["QualityTag"].append(line[10])
    return summary

def add_peaks(summary):
    # Extract peaks metrics and add them to summary
    summary["Peaks_number"] = []
    for peaks_file in snakemake.input.peaks:
        with open(peaks_file, "r") as f:
            summary["Peaks_number"].append(len(f.readlines()))
    return summary

if __name__ == "__main__":
    # Creation and saving of the dataframe
    summary = {"Sample": snakemake.params.samples}
    summary = add_mapping_reads(summary)
    summary = add_phantompeak(summary)
    summary = add_peaks(summary)
    summary_df = pd.DataFrame(summary)
    summary_df.to_csv(snakemake.output.path, sep="\t", index=False)