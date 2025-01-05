# Rules for read mapping and post-processing

# Build BWA-MEM2 index if needed
rule build_bwa_index:
    input:
        config["reference"]["fasta"]
    output:
        multiext(
            config["reference"]["index"],
            ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"
        )
    log:
        "logs/bwa_index/build.log"
    threads: 8
    shell:
        "bwa-mem2 index -p {config[reference][index]} {input} 2> {log}"

# Build Minimap2 index if needed
rule build_minimap2_index:
    input:
        config["reference"]["fasta"]
    output:
        config["reference"]["index"] + ".mmi"
    log:
        "logs/minimap2_index/build.log"
    threads: 8
    shell:
        "minimap2 -d {output} {input} 2> {log}"

# Map short reads with BWA-MEM2
rule map_short_reads:
    input:
        index = multiext(
            config["reference"]["index"],
            ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"
        ),
        r1 = lambda wildcards: config["samples"][wildcards.sample]["fq1"],
        r2 = lambda wildcards: config["samples"][wildcards.sample]["fq2"]
    output:
        temp("results/mapped/{sample}.raw.bam")
    params:
        extra = config["mapping"]["short_read"]["extra_params"],
        rg = lambda wildcards: f"@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tPL:illumina"
    log:
        "logs/bwa_mem2/{sample}.log"
    threads: config["threads"]
    shell:
        """
        bwa-mem2 mem -t {threads} {params.extra} -R '{params.rg}' \
            {config[reference][index]} {input.r1} {input.r2} 2> {log} | \
        samtools sort -@ {threads} -m 2G -O BAM -o {output}
        """

# Map long reads with Minimap2
rule map_long_reads:
    input:
        index = config["reference"]["index"] + ".mmi",
        reads = lambda wildcards: (
            config["samples"][wildcards.sample]["fq1"]
            if config["samples"][wildcards.sample]["type"] == "single"
            else config["samples"][wildcards.sample]["fq1"]
        )
    output:
        temp("results/mapped/{sample}.raw.bam")
    params:
        preset = lambda wildcards: (
            config["mapping"]["long_read_ont"]["preset"]
            if config["samples"][wildcards.sample]["platform"] == "ont"
            else config["mapping"]["long_read_pacbio"]["preset"]
        ),
        extra = lambda wildcards: (
            config["mapping"]["long_read_ont"]["extra_params"]
            if config["samples"][wildcards.sample]["platform"] == "ont"
            else config["mapping"]["long_read_pacbio"]["extra_params"]
        ),
        rg = lambda wildcards: (
            f"@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\t"
            f"PL:{config['samples'][wildcards.sample]['platform']}"
        )
    log:
        "logs/minimap2/{sample}.log"
    threads: config["threads"]
    shell:
        """
        minimap2 -ax {params.preset} -t {threads} {params.extra} \
            -R '{params.rg}' {input.index} {input.reads} 2> {log} | \
        samtools sort -@ {threads} -m 2G -O BAM -o {output}
        """

# Post-process BAM files
rule process_bam:
    input:
        "results/mapped/{sample}.raw.bam"
    output:
        bam = "results/mapped/{sample}.bam",
        bai = "results/mapped/{sample}.bam.bai"
    params:
        min_mapq = config["bam_processing"]["min_mapping_quality"],
        filter_supplementary = config["bam_processing"]["filter_supplementary"]
    log:
        "logs/process_bam/{sample}.log"
    threads: config["threads"]
    shell:
        """
        samtools view -@ {threads} -b \
            -q {params.min_mapq} \
            {"-F 0x800" if params.filter_supplementary else ""} \
            {input} > {output.bam} 2> {log}
            
        samtools index -@ {threads} {output.bam} {output.bai} 2>> {log}
        """
