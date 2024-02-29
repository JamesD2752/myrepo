"""
@author: James Damaso, COMP 383 001
"""
# myrepo
# When retrieving the transcriptomes, I used the SRA Toolkit which would download the SRA file in FASTQ format using fastq-dump, split the files into paired end files using --split-files
# and put the files in its original format using --origfmt
# I will be developing Track 2 


import os
import argparse

def run_pipeline_project(input_files, output_dir):
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    os.chdir(output_dir)

    # Step 2: Bowtie2 mapping
    bowtie2_index = "NC_006273.2"  # HCMV genome accession
    log_file = "PipelineProject.log"

    with open(log_file, "a") as log:
        log.write("Bowtie2 Mapping Results:\n")
        for file in input_files:
            # Run Bowtie2 mapping
            os.system(f"bowtie2 -x {bowtie2_index} -1 {file} -2 {file[:-6]}2.fastq -S {file[:-6]}.sam")

            # Count reads before filtering
            with open(file, "r") as f:
                reads_before = sum(1 for line in f) / 4

            # Count reads after filtering
            with open(f"{file[:-6]}2.fastq", "r") as f:
                reads_after = sum(1 for line in f) / 4

            # Write results to log
            log.write(f"{file[:-6]} had {int(reads_before)} read pairs before Bowtie2 filtering "
                      f"and {int(reads_after)} read pairs after.\n")

    # Step 3: Assembly (not implemented in this script)
    # Step 4: BLAST (not implemented in this script)
    # Step 5: Generate final output files (not implemented in this script)

    # Print summary to console
    print("Pipeline completed successfully.")

parser = argparse.ArgumentParser(description="Pipeline to process transcriptome reads.")
parser.add_argument("--input", nargs="+", help="Input FASTQ files", required=True)
parser.add_argument("--output", help="Output directory name", required=True)
args = parser.parse_args()

run_pipeline_project(args.input, args.output)
