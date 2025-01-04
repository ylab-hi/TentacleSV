# Mapping rules for both short and long reads

# Short-read mapping with bwa-mem2
rule map_short_reads:
    input:
        r1 = lambda wildcards: config["samples"][wildcards.sample]["fq1"],
        r2 = lambda wildcards: config["samples"][wildcards.sample]["fq2"],
        ref = config["reference"]
    output:
        "results/mapped/{sample}.bam"
    shell:
        """
        bwa-mem2 mem -t {threads} {input.ref} {input.r1} {input.r2} | \
        samtools sort -o {output}
        samtools index {output}
        """

# Long-read mapping with minimap2
rule map_long_reads:
    input:
        reads = lambda wildcards: config["samples"][wildcards.sample]["fq1"],
        ref = config["reference"]
    output:
        "results/mapped/{sample}.bam"
    shell:
        """
        minimap2 -ax map-ont -t {threads} {input.ref} {input.reads} | \
        samtools sort -o {output}
        samtools index {output}
        """
