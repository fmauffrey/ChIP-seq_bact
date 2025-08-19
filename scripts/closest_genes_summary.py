import pandas as pd
import numpy as np
import re

def load_peaks_table(file_path):
    """
    Load the peaks table.
    Keep only lines and columns with relevant information.
    Return a dataframe.
    """
    peaks_table = open(file_path, "r").readlines()
    peaks_table = [x for x in peaks_table if not x.startswith(("#", "\n", "chr\t"))]

    name = [x.split("\t")[9].strip("\n") for x in peaks_table]
    peak_location = [int(x.split("\t")[4]) for x in peaks_table]
    pileup = [float(x.split("\t")[5]) for x in peaks_table]
    fold_enrichment = [float(x.split("\t")[7]) for x in peaks_table]
    log10_qvalue = [float(x.split("\t")[8]) for x in peaks_table]

    extracted_data = {
        "Name": name,
        "Peak location": peak_location,
        "Pileup": pileup,
        "Fold enrichment": fold_enrichment,
        "Log10 Q-value": log10_qvalue
    }

    df = pd.DataFrame(extracted_data)

    return df

def load_annotation_table(file_path):
    """
    Load the annotation table.
    Keep only lines and columns with relevant information.
    """
    annotation_table = open(file_path, "r").readlines()
    annotation_table = [x for x in annotation_table if "\tgene\t" in x]

    names = [x.split("\t")[3] for x in annotation_table]
    gene_start = [int(x.split("\t")[8]) for x in annotation_table]
    gene_stop = [int(x.split("\t")[9]) for x in annotation_table]
    strand = [x.split("\t")[11] for x in annotation_table]
    info = [x.split("\t")[13] for x in annotation_table]
    products = []

    # Extract gene names
    gene_names = []
    for line in info:
        target_name = "NA"  # Default value
        for x in line.split(";"):
            if x.startswith("Name="):
                target_name = x.split("=")[1]
                break
        gene_names.append(target_name)
    
    # Extract product names if possible
    annotation_table = open(file_path, "r").read()
    for name in gene_names:
        pattern = re.compile(rf"\n.*\tCDS\t.*{name}.*\n", re.IGNORECASE)
        match = pattern.search(annotation_table)
        if match is not None:
            product_name = re.search(r"(?<=product=).*?(?=;)", match.group())
            products.append(product_name.group())
        else:
            products.append("NA")

    extracted_data = {
        "Name": names,
        "Gene name": gene_names,
        "Products": products,
        "Gene start": gene_start,
        "Gene stop": gene_stop,
        "Strand": strand,
    }

    df = pd.DataFrame(extracted_data)

    return df

def add_peak_distance(table):
    """
    Calculate the distance between a peak and a gene.
    Return 0 if the peak is located on the gene
    """
    table["Peak distance"] = np.where(
        table["Peak location"] < table["Gene start"],
        table["Gene start"] - table["Peak location"],
        np.where(
            table["Peak location"] > table["Gene stop"],
            table["Peak location"] - table["Gene stop"],
            0  # Peak is within the gene range
        )
    )

    return table

def add_TSS_state(table):
    """
    Check if the peak location could be considered as a TSS (located upstream).
    """
    table["Potential TSS"] = np.where(
        (table["Peak location"] < table["Gene start"]) & (table["Strand"] == "+"),
        "Yes",
        np.where(
            (table["Peak location"] > table["Gene stop"]) & (table["Strand"] == "-"),
            "Yes",
            "No"
        )
    )

    return table

if __name__ == "__main__":
    # Load peaks table and filter useless information
    peaks_table = load_peaks_table(snakemake.input.peaks_xls)
    annotation_table = load_annotation_table(snakemake.input.closest_genes)
    merged_table = pd.merge(peaks_table, annotation_table, on="Name", how="inner")
    merged_table = add_peak_distance(merged_table)
    merged_table = add_TSS_state(merged_table)

    # Write table
    merged_table.to_csv(snakemake.output.path, index=False, sep="\t")