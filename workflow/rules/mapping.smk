# Rules for read mapping and post-processing

# Build minimap2 index only if not provided
rule build_minimap2_index:
    input:
        ref = config["reference"]["fasta"]
    output:
        index = "results/genome_index/genome.mmi"
    log:
        "logs/minimap2_index/build.log"
    threads: config["threads"]
    shell:
        "minimap2 -d {output.index} {input.ref} 2> {log}"

# Get minimap2 index path (either user-provided or generated)
def get_minimap2_index(wildcards):
    if "minimap2_index" in config["reference"]:
        return config["reference"]["minimap2_index"]
    return "results/genome_index/genome.mmi"

# Map long reads (PacBio/ONT) with minimap2
rule map_long_reads:
    input:
        reads = lambda wildcards: config["samples"][wildcards.sample]["fq1"],
        index = get_minimap2_index
    output:
        temp("results/mapped/{sample}.raw.bam")
    params:
        # Set mapping mode based on platform
        preset = lambda wildcards: "map-hifi" if config["samples"][wildcards.sample]["platform"] == "pacbio" else "map-ont",
        rg = lambda wildcards: (
            f"@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\t"
            f"PL:{config['samples'][wildcards.sample]['platform']}"
        )
    log:
        "logs/minimap2/{sample}.log"
    threads: config["threads"]
    shell:
        """
        minimap2 -ax {params.preset} \
            --MD \
            -t {threads} \
            -Y \
            -R '{params.rg}' \
            {input.index} \
            {input.reads} 2> {log} | \
        samtools sort -@ {threads} -m 2G -O BAM -o {output}
        """

# Map short reads with BWA-MEM2 if needed
rule map_short_reads:
    input:
        index = multiext(
            config["reference"].get("bwa_index", "results/genome_index/genome"),
            ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"
        ),
        r1 = lambda wildcards: config["samples"][wildcards.sample]["fq1"],
        r2 = lambda wildcards: config["samples"][wildcards.sample]["fq2"]
    output:
        temp("results/mapped/{sample}.raw.bam")
    params:
        rg = lambda wildcards: f"@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tPL:illumina"
    log:
        "logs/bwa_mem2/{sample}.log"
    threads: config["threads"]
    shell:
        """
        bwa-mem2 mem \
            -t {threads} \
            -M \
            -Y \
            -R '{params.rg}' \
            {input.index} \
            {input.r1} \
            {input.r2} 2> {log} | \
        samtools sort -@ {threads} -m 2G -O BAM -o {output}
        """

# Process and index BAM files
rule process_bam:
    input:
        "results/mapped/{sample}.raw.bam"
    output:
        bam = "results/mapped/{sample}.bam",
        bai = "results/mapped/{sample}.bam.bai"
    log:
        "logs/process_bam/{sample}.log"
    params:
        min_mapq = 20,  # Default mapping quality threshold
        extra_flags = "-F 0x800"  # Filter supplementary alignments
    threads: config["threads"]
    shell:
        """
        # Filter BAM file
        samtools view \
            -@ {threads} \
            -b \
            -q {params.min_mapq} \
            {params.extra_flags} \
            {input} > {output.bam} 2> {log}
        
        # Index filtered BAM
        samtools index -@ {threads} {output.bam} {output.bai} 2>> {log}
        """
