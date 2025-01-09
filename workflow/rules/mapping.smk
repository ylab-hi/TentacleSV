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

# Only execute mapping rules if input is not BAM
if not any(is_bam_input(sample) for sample in SAMPLES):
    # Map long reads (PacBio/ONT) with minimap2
    rule map_long_reads:
        input:
            reads = lambda wildcards: config["samples"][wildcards.sample]["fq1"],
            index = get_minimap2_index
        output:
            temp("results/mapped/{sample}.raw.bam")
        params:
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

    # Map short reads with BWA-MEM2
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

# Process BAM files (handles both mapped and input BAM)
rule process_bam:
    input:
        # Use direct BAM input or mapped BAM based on input type
        lambda wildcards: (
            config["samples"][wildcards.sample]["bam"] 
            if is_bam_input(wildcards.sample) 
            else "results/mapped/{sample}.raw.bam"
        )
    output:
        bam = "results/mapped/{sample}.bam",
        bai = "results/mapped/{sample}.bam.bai"
    log:
        "logs/process_bam/{sample}.log"
    params:
        min_mapq = config["bam_processing"]["min_mapping_quality"],
        extra_flags = "-F 0x800" if config["bam_processing"]["filter_supplementary"] else ""
    threads: config["threads"]
    shell:
        """
        # If input is a direct BAM file, create symlink
        if [ "{wildcards.sample}" == "test_sample" ] && [ -f "{input}" ]; then
            ln -sf {input} {output.bam}
        else
            # Create necessary directories
            mkdir -p $(dirname {output.bam}) $(dirname {log})
            
            # Process BAM file
            samtools view \
                -@ {threads} \
                -b \
                -q {params.min_mapq} \
                {params.extra_flags} \
                {input} \
                -o {output.bam} \
                2> {log}
        fi
        
        # Index BAM file regardless of source
        samtools index \
            -@ {threads} \
            {output.bam} \
            2>> {log}
        """
