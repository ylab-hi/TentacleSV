# Common functions and configurations for TentacleSV workflow
import os
from snakemake.utils import validate
from typing import Dict, List, Union

# Get all sample IDs
SAMPLES = list(config["samples"].keys())

# Determine sequencing type and appropriate callers
SEQ_TYPE = config["seq_type"]
if SEQ_TYPE == "short":
    CALLERS = config["callers"]["short_read"]
elif SEQ_TYPE == "long_ont":
    CALLERS = config["callers"]["long_read_ont"]
elif SEQ_TYPE == "long_pacbio":
    CALLERS = config["callers"]["long_read_pacbio"]
else:
    raise ValueError(f"Unknown sequencing type: {SEQ_TYPE}")

def get_input_files(wildcards: snakemake.io.Wildcards) -> Union[str, List[str]]:
    """
    Get input files for a given sample based on sample configuration.
    
    Args:
        wildcards: Snakemake wildcards containing sample information
        
    Returns:
        str or List[str]: Path(s) to input file(s)
    """
    sample_config = config["samples"][wildcards.sample]
    if sample_config["type"] == "bam":
        return sample_config["bam"]
    elif sample_config["type"] == "paired":
        return [sample_config["fq1"], sample_config["fq2"]]
    elif sample_config["type"] == "single":
        return sample_config["fq1"]
    else:
        raise ValueError(f"Unknown sample type for {wildcards.sample}")

def get_mapping_params(wildcards: snakemake.io.Wildcards) -> Dict[str, str]:
    """
    Get mapping parameters based on sequencing type and platform.
    
    Args:
        wildcards: Snakemake wildcards containing sample information
        
    Returns:
        Dict[str, str]: Mapping tool and parameters
    """
    sample_config = config["samples"][wildcards.sample]
    platform = sample_config["platform"]
    
    if platform == "illumina":
        params = config["mapping"]["short_read"]
    elif platform == "ont":
        params = config["mapping"]["long_read_ont"]
    elif platform == "pacbio":
        params = config["mapping"]["long_read_pacbio"]
    else:
        raise ValueError(f"Unknown platform: {platform}")
    
    return params

def get_final_output() -> List[str]:
    """
    Generate list of final output files expected from the workflow.
    
    Returns:
        List[str]: Paths to expected output files
    """
    outputs = []
    
    # Add merged VCF for each sample
    outputs.extend(expand(
        "results/merged/{sample}.vcf",
        sample=SAMPLES
    ))
    
    # Add UpSet plot if enabled
    if config["output"]["upset_plot"]:
        outputs.append("results/plots/upset.png")
    
    return outputs

# Common rules that might be needed by multiple modules
rule index_bam:
    input:
        "results/mapped/{sample}.bam"
    output:
        "results/mapped/{sample}.bam.bai"
    log:
        "logs/samtools/{sample}.index.log"
    threads: 4
    shell:
        "samtools index -@ {threads} {input} {output} 2> {log}"

rule remove_duplicates:
    input:
        bam = "results/mapped/{sample}.bam",
        bai = "results/mapped/{sample}.bam.bai"
    output:
        bam = "results/mapped/{sample}.rmdup.bam",
        bai = "results/mapped/{sample}.rmdup.bam.bai"
    log:
        "logs/samtools/{sample}.rmdup.log"
    threads: 4
    shell:
        """
        samtools rmdup {input.bam} {output.bam} 2> {log}
        samtools index -@ {threads} {output.bam} {output.bai} 2>> {log}
        """
