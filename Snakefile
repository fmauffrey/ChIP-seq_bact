rule all:
    input:
        "data/processed_data.txt"

rule process_data:
    input:
        "data/raw_data.txt"
    output:
        "data/processed_data.txt"
    script:
        "scripts/process_data.py"