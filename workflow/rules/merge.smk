# Rules for OctopusV correction and merging

# Convert VCF to SVCF and correct BND events
rule correct_vcf:
    input:
        vcf = "results/sv_calls/{caller}/{sample}.vcf"
    output:
        svcf = "results/corrected/{caller}/{sample}.svcf"
    params:
        pos_tolerance = config["octopusv"]["correct"]["pos_tolerance"]
    log:
        "logs/octopusv/correct_{caller}_{sample}.log"
    shell:
        """
        octopusv correct \
            -i {input.vcf} \
            -o {output.svcf} \
            --pos-tolerance {params.pos_tolerance} \
            2> {log}
        """

# Merge corrected SVCFs from all callers for a sample
rule merge_sample_calls:
    input:
        svcfs = lambda wildcards: expand(
            "results/corrected/{caller}/{sample}.svcf",
            caller=CALLERS,
            sample=wildcards.sample
        )
    output:
        merged_svcf = "results/merged/{sample}.svcf"
    params:
        min_support = config["octopusv"]["merge"]["min_support"],
        max_distance = config["octopusv"]["merge"]["max_distance"],
        max_length_ratio = config["octopusv"]["merge"]["max_length_ratio"],
        min_jaccard = config["octopusv"]["merge"]["min_jaccard"],
        tra_delta = config["octopusv"]["merge"]["tra_delta"],
        tra_min_overlap = config["octopusv"]["merge"]["tra_min_overlap"],
        tra_strand_consistency = "--tra-strand-consistency" if config["octopusv"]["merge"]["tra_strand_consistency"] else ""
    log:
        "logs/octopusv/merge_{sample}.log"
    shell:
        """
        octopusv merge \
            -i {input.svcfs} \
            -o {output.merged_svcf} \
            --min-support {params.min_support} \
            --max-distance {params.max_distance} \
            --max-length-ratio {params.max_length_ratio} \
            --min-jaccard {params.min_jaccard} \
            --tra-delta {params.tra_delta} \
            --tra-min-overlap {params.tra_min_overlap} \
            {params.tra_strand_consistency} \
            2> {log}
        """

# Convert merged SVCF back to VCF
rule svcf_to_vcf:
    input:
        svcf = "results/merged/{sample}.svcf"
    output:
        vcf = "results/merged/{sample}.vcf"
    log:
        "logs/octopusv/svcf2vcf_{sample}.log"
    shell:
        """
        octopusv svcf2vcf \
            -i {input.svcf} \
            -o {output.vcf} \
            2> {log}
        """

# Generate UpSet plot if enabled
rule generate_upset_plot:
    input:
        svcfs = lambda wildcards: expand(
            "results/corrected/{caller}/{sample}.svcf",
            caller=CALLERS,
            sample=wildcards.sample
        )
    output:
        plot = "results/plots/{sample}_upset.png"
    params:
        merged_svcf = "results/merged/{sample}.svcf"
    log:
        "logs/octopusv/upset_plot_{sample}.log"
    shell:
        """
        octopusv merge \
            -i {input.svcfs} \
            -o {params.merged_svcf} \
            --upsetr \
            --upsetr-output {output.plot} \
            2> {log}
        """
