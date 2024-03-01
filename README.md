"""
@author: James Damaso, COMP 383 001
"""
# When retrieving the transcriptomes, I used the SRA Toolkit which would download the SRA file in FASTQ format using fastq-dump, split the files into paired-end files using --split-files, and put the files in its original format using --origfmt.
# I will be developing Track 1 

Required Software
Linux/Unix
Python3
Biopython
Kallisto
R

Entire code below

used this in the command line to start algorithm:
python PipelineProject_Damaso.py NC_006273.2 SRR5660030_1 SRR5660030_2 SRR5660033_1 SRR5660033_2 SRR5660044_1 SRR5660044_2 SRR5660045_1 SRR5660045_2 --input input.gb --output-dir PipelineProject_James_Damaso

import os
import subprocess
import argparse
import statistics
from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML

# Function to remove null bytes from a file
def remove_null_bytes(input_file, output_file):
    """
    Removes null bytes from a file and writes the cleaned content to another file.
    Parameters:
        input_file (str): Path to the input file.
        output_file (str): Path to the output file.
    """
    with open(input_file, 'rb') as f:
        data = f.read().replace(b'\x00', b'')
    with open(output_file, 'wb') as f:
        f.write(data)

# Function to extract CDS (coding sequences) from a GenBank file and save them in a FASTA file
def extract_cds(input_gb_file, output_fasta_file):
    """
    Extracts CDS features from a GenBank file and writes them to a FASTA file.
    Parameters:
        input_gb_file (str): Path to the input GenBank file.
        output_fasta_file (str): Path to the output FASTA file.
    Returns:
        int: Number of CDS features extracted.
    """
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

# Function to run BLAST search using NCBI's online service
def run_blast(input_fasta, output_xml):
    """
    Runs BLAST search using NCBI's online service.
    Parameters:
        input_fasta (str): Path to the input FASTA file.
        output_xml (str): Path to the output XML file containing BLAST results.
    """
    blast_db = "path_to_betaherpesvirinae_db"  # Replace with the actual path to your local BLAST database
    result_handle = NCBIWWW.qblast("blastx", blast_db, input_fasta)
    with open(output_xml, "w") as out_file:
        out_file.write(result_handle.read())
    result_handle.close()

# Function to parse BLAST results and write them into a log file
def parse_blast_results(xml_file, log_file):
    """
    Parses BLAST results from an XML file and writes them into a log file.
    Parameters:
        xml_file (str): Path to the input XML file containing BLAST results.
        log_file (str): Path to the output log file.
    """
    with open(log_file, "a") as log:
        log.write("target_id\ttest_stat\tpval\tqval\n")
        result_handle = open(xml_file)
        blast_records = NCBIXML.parse(result_handle)
        for blast_record in blast_records:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    log.write(f"{alignment.accession}\t{hsp.identities / hsp.align_length * 100:.2f}\t{hsp.align_length}\t{hsp.query_start}\t{hsp.query_end}\t{hsp.sbjct_start}\t{hsp.sbjct_end}\t{hsp.bits}\t{hsp.expect}\t{alignment.title}\n")
        result_handle.close()

# Function to build a kallisto index for RNA-seq quantification
def build_kallisto_index(input_fasta_file, output_directory):
    """
    Builds a kallisto index for RNA-seq quantification.
    Parameters:
        input_fasta_file (str): Path to the input FASTA file.
        output_directory (str): Path to the output directory for storing the index file.
    Returns:
        str: Path to the generated kallisto index file.
    """
    index_file = os.path.join(output_directory, "output_index.idx")
    subprocess.run(["kallisto", "index", "-i", index_file, input_fasta_file], check=True)
    return index_file

# Function to quantify transcript abundance using kallisto
def quantify_tpm(sample, index_file, output_directory):
    """
    Quantifies transcript abundance using kallisto.
    Parameters:
        sample (str): Sample name.
        index_file (str): Path to the kallisto index file.
        output_directory (str): Path to the output directory for storing the quantification results.
    """
    fastq_file = f"{sample}.fastq"
    output_dir = os.path.join(output_directory, sample)
    os.makedirs(output_dir, exist_ok=True)  # Create output directory
    subprocess.run(["kallisto", "quant", "-i", index_file, "-o", output_dir, "--single", "-l", "200", "-s", "20", fastq_file], check=True)

# Function to extract condition information from a sample file
def extract_condition_from_file(sample_file):
    """
    Extracts condition information from a sample file.
    Parameters:
        sample_file (str): Path to the sample file.
    Returns:
        str: Condition information extracted from the sample file.
    """
    condition = None
    with open(sample_file, "r") as f:
        for line in f:
            if line.startswith("@"):  # FASTQ header starts with '@'
                parts = line.strip().split()  # Split header line by whitespace
                for part in parts:
                    if part.startswith("#Condition:"):  # Assuming condition is stored as '#Condition:'
                        condition = part.split(":")[1].strip()
                        break
                break  # Breaks once the header line has been processed
    return condition

# Function to calculate statistics from kallisto output
def calculate_tpm_stats(sample_dir, condition):
    """
    Calculates statistics from kallisto output.
    Parameters:
        sample_dir (str): Path to the directory containing kallisto output for a sample.
        condition (str): Condition information associated with the sample.
    Returns:
        tuple: Minimum TPM, median TPM, mean TPM, and maximum TPM.
    """
    tpm_values = []
    with open(os.path.join(sample_dir, "abundance.tsv"), "r") as f:
        next(f)  # skip header
        for line in f:
            tpm_values.append(float(line.split()[4]))
    min_tpm = min(tpm_values)
    med_tpm = statistics.median(tpm_values)
    mean_tpm = statistics.mean(tpm_values)
    max_tpm = max(tpm_values)
    return min_tpm, med_tpm, mean_tpm, max_tpm

# Function to run sleuth R script for differential expression analysis
def run_sleuth(output_dir, sample_dir):
    """
    Runs sleuth R script for differential expression analysis.
    Parameters:
        output_dir (str): Path to the output directory.
        sample_dir (str): Path to the directory containing kallisto output for a sample.
    """
    script_path = os.path.join(os.path.dirname(__file__), "sleuth_script.R")
    try:
        subprocess.run(["Rscript", script_path, sample_dir], check=True)
    except subprocess.CalledProcessError as e:
        print("Error running R script:", e)
        return
    sleuth_results_file = os.path.join(sample_dir, "topten.txt")
    if os.path.exists(sleuth_results_file):
        print("Sleuth results file found:", sleuth_results_file)
        with open(sleuth_results_file, "r") as f:
            print("Target ID\tTest Stat\tP-value\tQ-value")
            for line in f:
                target_id, test_stat, pval, qval = line.strip().split("\t")
                print(target_id, test_stat, pval, qval)
    else:
        print("No Sleuth results found for", sample_dir)

# Main function that drives the entire pipeline
def main(args):
    input_gb_file = args.input
    output_fasta_file = os.path.join(args.output_dir, "output.fasta")
    log_file = os.path.join(args.output_dir, "PipelineProject.log")
    samples = args.samples

    os.makedirs(args.output_dir, exist_ok=True)

    num_cds = extract_cds(input_gb_file, output_fasta_file)

    index_file = build_kallisto_index(output_fasta_file, args.output_dir)

    with open(log_file, "w") as log:
        log.write(f"The HCMV genome ({args.accession}) has {num_cds} CDS.\n\n")
        log.write("sample\tcondition\tmin_tpm\tmed_tpm\tmean_tpm\tmax_tpm\n")
        for sample in samples:
            condition = extract_condition_from_file(f"{sample}.fastq")  # Extracts condition from sample file
            quantify_tpm(sample, index_file, args.output_dir)  # Quantify TPM
            min_tpm, med_tpm, mean_tpm, max_tpm = calculate_tpm_stats(os.path.join(args.output_dir, sample), condition)
            log.write(f"{sample}\t{condition}\t{min_tpm}\t{med_tpm}\t{mean_tpm}\t{max_tpm}\n")

        # Writes sleuth header to the log file
        log.write("target_id\ttest_stat\tpval\tqval\n")

    # Runs BLAST search
    blast_output_xml = os.path.join(args.output_dir, "blast_output.xml")
    run_blast(output_fasta_file, blast_output_xml)

    # Parses BLAST results and append to log file
    parse_blast_results(blast_output_xml, log_file)

    # Runs sleuth to find differentially expressed genes
    run_sleuth(args.output_dir, sample_dir=os.path.join(args.output_dir, samples[0]))

# Parsing command-line arguments
parser = argparse.ArgumentParser(description="Build kallisto index for a given NCBI accession number and quantify TPM for each sample.")
parser.add_argument("accession", type=str, help="NCBI accession number")
parser.add_argument("--input", type=str, default="input.gb", help="Input GenBank file (default: input.gb)")
parser.add_argument("--output-dir", type=str, default="PipelineProject_James_Damaso", help="Output directory (default: PipelineProject_James_Damaso)")
parser.add_argument("samples", nargs="+", help="List of sample names (SRR numbers)")
args = parser.parse_args()

# Calling the main function with parsed arguments
main(args)
