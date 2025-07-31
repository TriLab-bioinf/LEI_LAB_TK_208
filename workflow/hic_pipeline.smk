# rRNA-depletion pipeline based on SortMeRNA
# Contact: Hernan Lorenzi
# Contact email: hernan.lorenzi@nih.gov
# Workflow asumes reads are already trimmed and deduplicated (based on UMIs) 
# Workflow requires to configure the config.yml file to set folder containing reads and rRNA fasta file and database folder.
# Workflow also requires to to edit samplesheet.csv file with sample information 

import os
import glob
import pandas as pd

# Read config file
configfile: "./config/config.yaml"
adapters: str = config["adapters"]
reads_dir: str = config["reads"]

# Read sample data from samplesheet and skip comments
metadata = pd.read_csv(config["metadata"], comment='#', sep=',', header=0, dtype=str)

# Create dictionaries from metadata for @RG line
# metadata is imported from ./config/samplesheet.csv file
fq_1: dict = {s:fq1 for s, fq1 in zip(metadata['sample_ID'], metadata['fastq_1'])}
fq_2: dict = {s:fq2 for s, fq2 in zip(metadata['sample_ID'], metadata['fastq_2'])}
samples: list = list(metadata['sample_ID'])

# Set what rules to run locally
localrules: all 

rule all:
    # IMPORTANT: output file for all rules has to match the name specified in the output file
    # and not include suffixes that the command use might add to it.
    input: "results/3-multiqc/multiqc_report.html"  

# 1- RTrim reads with fastp
if config["trimming"] == True:
    trim_output = "results/1-trimming"
    # 1- Trim reads with fastp
    include: "rules/trim_reads.smk"
else:
    trim_output = "data/reads"

# 6- Run multiqc report
include: "rules/multiqc.smk"
