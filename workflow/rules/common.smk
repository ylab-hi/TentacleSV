# Common functions and variables for TentacleSV workflow
import os
from snakemake.utils import validate
from snakemake.utils import min_version

# Set minimum snakemake version
min_version("6.0.0")

# Create required directories
for d in ["data", "results/mapped", "results/sv_calls", "results/corrected", 
          "results/merged", "logs/mapping", "logs/sv_calls"]:
    os.makedirs(d, exist_ok=True)

# Get all sample IDs
SAMPLES = list(config["samples"].keys())

# Helper function to get input files
def get_input_files(wildcards):
    sample_config = config["samples"][wildcards.sample]
    if sample_config["type"] == "single":
        return sample_config["fq1"]
    elif sample_config["type"] == "paired":
        return [sample_config["fq1"], sample_config["fq2"]]
    elif sample_config["type"] == "bam":
        return sample_config["bam"]
    else:
        raise ValueError(f"Unknown sample type for {wildcards.sample}")

# Get appropriate callers based on sequencing type
def get_callers():
    seq_type = config["seq_type"]
    if seq_type == "short":
        return config["callers"].get("short_read", [])
    elif seq_type == "long_pacbio":
        return config["callers"].get("long_read_pacbio", [])
    elif seq_type == "long_ont":
        return config["callers"].get("long_read_ont", [])
    else:
        raise ValueError(f"Unknown sequencing type: {seq_type}")

# Set active callers
CALLERS = get_callers()
