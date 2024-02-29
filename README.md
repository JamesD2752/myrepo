"""
@author: James Damaso, COMP 383 001
"""
# myrepo
# When retrieving the transcriptomes, I used the SRA Toolkit which would download the SRA file in FASTQ format using fastq-dump, split the files into paired end files using --split-files
# and put the files in its original format using --origfmt
# I will be developing Track 1 

import os
import subprocess
import argparse
import statistics
from Bio import SeqIO

def extract_cds(input_gb_file, output_fasta_file):
    cds_records = []
    with open(input_gb_file, "r") as handle:
        for record in SeqIO.parse(handle, "genbank"):
            for feature in record.features:
                if feature.type == "CDS":
                    cds_records.append(feature)

    with open(output_fasta_file, "w") as output_handle:
        for cds in cds_records:
            protein_id = cds.qualifiers.get('protein_id', [''])[0]
            if protein_id:
                output_handle.write(f">{protein_id}\n{cds.location.extract(record).seq}\n")

    return len(cds_records)

def build_kallisto_index(input_fasta_file):
    index_file = os.path.splitext(input_fasta_file)[0] + "_index.idx"
    subprocess.run(["kallisto", "index", "-i", index_file, input_fasta_file], check=True, stdout=subprocess.PIPE)
    return index_file

def quantify_tpm(sample, index_file):
    subprocess.run(["kallisto", "quant", "-i", index_file, "-o", sample, f"{sample}_reads.fastq"], check=True, stdout=subprocess.PIPE)

def calculate_tpm_stats(sample_dir):
    tpm_values = []
    with open(os.path.join(sample_dir, "abundance.tsv"), "r") as f:
        next(f)  # skip header
        for line in f:
            tpm_values.append(float(line.split()[4]))
    return min(tpm_values), statistics.median(tpm_values), statistics.mean(tpm_values), max(tpm_values)

def main(args):
    input_gb_file = args.input
    output_fasta_file = args.output
    log_file = args.log
    samples = args.samples

    # Extract CDS and write to FASTA file
    num_cds = extract_cds(input_gb_file, output_fasta_file)

    # Build kallisto index
    index_file = build_kallisto_index(output_fasta_file)

    # Quantify TPM for each sample
    for sample in samples:
        quantify_tpm(sample, index_file)

    # Write log information
    with open(log_file, "w") as log:
        log.write("sample\tcondition\tmin_tpm\tmed_tpm\tmean_tpm\tmax_tpm\n")
        for sample in samples:
            min_tpm, med_tpm, mean_tpm, max_tpm = calculate_tpm_stats(os.path.join(sample, "abundance.tsv"))
            log.write(f"{sample}\t{args.condition}\t{min_tpm}\t{med_tpm}\t{mean_tpm}\t{max_tpm}\n")

parser = argparse.ArgumentParser(description="Build kallisto index for a given NCBI accession number and quantify TPM for each sample.")
parser.add_argument("accession", type=str, help="NCBI accession number")
parser.add_argument("--input", type=str, default="input.gb", help="Input GenBank file (default: input.gb)")
parser.add_argument("--output", type=str, default="output.fasta", help="Output FASTA file (default: output.fasta)")
parser.add_argument("--log", type=str, default="log.txt", help="Log file (default: log.txt)")
parser.add_argument("--condition", type=str, default="unknown", help="Condition of the samples (default: unknown)")
parser.add_argument("samples", nargs="+", help="List of sample names (SRR numbers)")
args = parser.parse_args()
main(args)

