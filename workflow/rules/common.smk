# Common functions and variables for TentacleSV workflow
import os

# Helper function to get input files based on sequencing type
def get_input_files(wildcards):
    """Return input files for a given sample"""
    if config["samples"][wildcards.sample].get("bam"):
        return config["samples"][wildcards.sample]["bam"]
    else:
        return [
            config["samples"][wildcards.sample]["fq1"],
            config["samples"][wildcards.sample]["fq2"]
        ]

# Get all sample IDs
SAMPLES = list(config["samples"].keys())

# Select SV callers based on sequencing type
CALLERS = config["callers"]["short_read"] if config["seq_type"] == "short" else config["callers"]["long_read"]
