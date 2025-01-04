# Rules for OctopusV correction and merging

# Correct BND events in each VCF
rule correct_vcf:
    input:
        vcf = "results/sv_calls/{caller}/{sample}.vcf"
    output:
        svcf = "results/corrected/{caller}/{sample}.svcf"
    params:
        pos_tolerance = config["correct"]["pos_tolerance"]
    shell:
        """
        octopusv correct \
            -i {input.vcf} \
            -o {output.svcf} \
            --pos-tolerance {params.pos_tolerance}
        """

# Merge all corrected SVCFs for a sample
rule merge_sample_calls:
    input:
        svcfs = expand("results/corrected/{caller}/{{sample}}.svcf", caller=CALLERS)
    output:
        merged = "results/merged/{sample}.vcf"
    params:
        min_support = config["merge"]["min_support"],
        max_distance = config["merge"]["max_distance"],
        max_length_ratio = config["merge"]["max_length_ratio"],
        min_jaccard = config["merge"]["min_jaccard"]
    shell:
        """
        octopusv merge \
            -i {input.svcfs} \
            -o {output.merged} \
            --min-support {params.min_support} \
            --max-distance {params.max_distance} \
            --max-length-ratio {params.max_length_ratio} \
            --min-jaccard {params.min_jaccard}
        """
