"""
Rules for structural variant calling.
This module contains rules for running different SV callers on mapped reads.
"""

# Helper function to get tandem repeats file path
def get_tandem_repeats():
    """
    Get tandem repeats file path if available.
    Returns:
        str: Path to tandem repeats file or empty string
    """
    return config["reference"].get("tandem_repeats", "")

# Helper function to construct tool-specific arguments
def get_tr_argument(wildcards, tool):
    """
    Get tool-specific tandem repeats argument.
    Args:
        wildcards: Snakemake wildcards
        tool: Name of the tool
    Returns:
        str: Command line argument for tandem repeats
    """
    tr_path = get_tandem_repeats()
    if not tr_path:
        return ""
    
    if tool == "sniffles":
        return f"--tandem-repeats {tr_path}"
    elif tool == "pbsv":
        return f"--tandem-repeats {tr_path}"
    return ""

# Short-read callers
rule run_manta:
    input:
        bam = "results/mapped/{sample}.bam",
        bai = "results/mapped/{sample}.bam.bai",
        ref = config["reference"]["fasta"]
    output:
        vcf = "results/sv_calls/manta/{sample}.vcf"
    container:
        config["containers"]["manta"]
    log:
        "logs/manta/{sample}.log"
    threads: config["threads"]
    shell:
        """
        # Create necessary directories
        mkdir -p logs/manta results/sv_calls/manta
        TMPDIR=$(mktemp -d)
        
        # Configure Manta
        configManta.py \
            --bam {input.bam} \
            --referenceFasta {input.ref} \
            --runDir $TMPDIR \
            2> {log}
        
        # Run workflow
        $TMPDIR/runWorkflow.py -j {threads} 2>> {log}
        
        # Move and process results
        mv $TMPDIR/results/variants/diploidSV.vcf.gz {output.vcf}.gz && \
        gunzip {output.vcf}.gz && \
        rm -rf $TMPDIR
        """


rule run_delly:
    """
    Run Delly structural variant caller.
    Specializes in deletion, duplication, inversion and translocation detection.
    """
    input:
        bam = "results/mapped/{sample}.bam",
        bai = "results/mapped/{sample}.bam.bai",
        ref = config["reference"]["fasta"]
    output:
        vcf = "results/sv_calls/delly/{sample}.vcf"
    params:
        min_mapq = config["delly"]["min_mapping_quality"],
        pe_support = config["delly"]["min_paired_end_support"],
        max_insert = config["delly"]["max_insert_size"],
        min_size = config["delly"]["min_size"]
    conda:
        "../../envs/environment.full.yaml"
    log:
        "logs/delly/{sample}.log"
    threads: config["threads"]
    shell:
        """
        # Run DELLY with configured parameters
        delly call \
            -g {input.ref} \
            -q {params.min_mapq} \
            -s {params.pe_support} \
            -n {params.max_insert} \
            -m {params.min_size} \
            -o {output.vcf}.tmp \
            {input.bam} \
            2> {log}
            
        # Filter for PASS variants
        bcftools filter -i 'FILTER="PASS"' \
            {output.vcf}.tmp > {output.vcf} \
            2>> {log}
            
        rm {output.vcf}.tmp
        """

rule run_lumpy:
    input:
        bam = "results/mapped/{sample}.bam",
        bai = "results/mapped/{sample}.bam.bai"
    output:
        vcf = "results/sv_calls/lumpy/{sample}.vcf"
    log:
        "logs/lumpy/{sample}.log"
    container:
        config["containers"]["lumpy"]
    threads: 
        config["threads"]
    shell:
        """
        # Create necessary directories
        mkdir -p results/sv_calls/lumpy logs/lumpy
        
        # Run LUMPY
        lumpyexpress \
            -B {input.bam} \
            -o {output.vcf} \
            2> {log}
        """


rule run_svaba:
    """
    Run SvABA structural variant caller
    """
    input:
        bam = "results/mapped/{sample}.bam",
        bai = "results/mapped/{sample}.bam.bai",
        ref = config["reference"]["fasta"]
    output:
        vcf = "results/sv_calls/svaba/{sample}.vcf"
    conda:
        "../../envs/environment.full.yaml"
    log:
        "logs/svaba/{sample}.log"
    threads: config["threads"]
    shell:
        """
        # Create output directories
        mkdir -p results/sv_calls/svaba logs/svaba
        
        svaba run \
            -t $(readlink -f {input.bam}) \
            -p {threads} \
            -a TentacleSV \
            -G $(readlink -f {input.ref}) \
            2> {log}
        
        # Move the SV VCF to the final location
        mv TentacleSV.svaba.sv.vcf {output.vcf}
        
        # Clean up intermediate files
        rm -f TentacleSV.*
        rm -rf results/sv_calls/svaba/work_*
        """


rule run_cutesv:
    """
    Run CuteSV structural variant caller.
    Optimized for long-read sequencing data.
    """
    input:
        bam = "results/mapped/{sample}.bam",
        bai = "results/mapped/{sample}.bam.bai",
        ref = config["reference"]["fasta"]
    output:
        vcf = "results/sv_calls/cutesv/{sample}.vcf"
    params:
        temp_dir = "results/sv_calls/cutesv/tmp_{sample}",
        min_support = config["cutesv"]["min_support"],
        min_size = config["cutesv"]["min_size"],
        max_size = config["cutesv"]["max_size"],
        min_mapq = config["cutesv"]["min_mapq"],
        signal_distance = config["cutesv"]["signal_distance"],
        diff_ratio_merging_ins = config["cutesv"]["diff_ratio_merging_ins"],
        diff_ratio_merging_del = config["cutesv"]["diff_ratio_merging_del"]
    conda:
        "../../envs/environment.full.yaml"
    log:
        "logs/cutesv/{sample}.log"
    threads: config["threads"]
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 8000
    shell:
        """
        # Create temporary directory
        rm -rf {params.temp_dir}
        mkdir -p {params.temp_dir}
        
        # Run CuteSV with configured parameters
        cuteSV \
            --threads {threads} \
            --genotype \
            --report_readid \
            --min_support {params.min_support} \
            --min_size {params.min_size} \
            --max_size {params.max_size} \
            --min_mapq {params.min_mapq} \
            --signal_distance {params.signal_distance} \
            --diff_ratio_merging_INS {params.diff_ratio_merging_ins} \
            --diff_ratio_merging_DEL {params.diff_ratio_merging_del} \
            {input.bam} \
            {input.ref} \
            {output.vcf} \
            {params.temp_dir} \
            2> {log}
        
        # Clean up
        rm -rf {params.temp_dir}
        
        # Verify output
        if [ ! -s {output.vcf} ]; then
            echo "CuteSV completed but output file is empty" >> {log}
            exit 1
        fi
        """

rule run_sniffles:
    """
    Run Sniffles structural variant caller.
    Specialized in detecting SVs from long-read sequencing data.
    """
    input:
        bam = "results/mapped/{sample}.bam",
        bai = "results/mapped/{sample}.bam.bai",
        ref = config["reference"]["fasta"]
    output:
        vcf = "results/sv_calls/sniffles/{sample}.vcf"
    params:
        tr_arg = lambda wildcards: get_tr_argument(wildcards, "sniffles"),
        min_support = config["sniffles"]["min_support"],
        min_length = config["sniffles"]["min_length"],
        max_num_splits = config["sniffles"]["max_num_splits"],
        max_distance = config["sniffles"]["max_distance"],
        min_alignment_length = config["sniffles"]["min_alignment_length"]
    conda:
        "../../envs/environment.full.yaml"
    log:
        "logs/sniffles/{sample}.log"
    threads: config["threads"]
    shell:
        """
        sniffles \
            --input {input.bam} \
            --vcf {output.vcf} \
            --reference {input.ref} \
            {params.tr_arg} \
            --minsupport {params.min_support} \
            --minsvlen {params.min_length} \
            --max-splits {params.max_num_splits} \
            --max-distance {params.max_distance} \
            --min-alignment-length {params.min_alignment_length} \
            --threads {threads} \
            2> {log}
        """

rule run_svim:
    """
    Run SVIM structural variant caller.
    Specialized for long-read data analysis.
    """
    input:
        bam = "results/mapped/{sample}.bam",
        bai = "results/mapped/{sample}.bam.bai",
        ref = config["reference"]["fasta"]
    output:
        vcf = "results/sv_calls/svim/{sample}.vcf"
    params:
        min_mapq = config["svim"]["min_mapq"],
        min_sv_size = config["svim"]["min_sv_size"],
        max_sv_size = config["svim"]["max_sv_size"],
        min_read_length = config["svim"]["min_read_length"],
        gap_tolerance = config["svim"]["segment_gap_tolerance"],
        overlap_tolerance = config["svim"]["segment_overlap_tolerance"]
    conda:
        "../../envs/environment.full.yaml"
    log:
        "logs/svim/{sample}.log"
    shell:
        """
        # Create temporary directory
        TMPDIR=$(mktemp -d)
        
        # Run SVIM with configured parameters
        svim alignment \
            --min_mapq {params.min_mapq} \
            --min_sv_size {params.min_sv_size} \
            --max_sv_size {params.max_sv_size} \
            --min_read_length {params.min_read_length} \
            --segment_gap_tolerance {params.gap_tolerance} \
            --segment_overlap_tolerance {params.overlap_tolerance} \
            $TMPDIR \
            {input.bam} \
            {input.ref} \
            2> {log}
            
        # Move results
        mv $TMPDIR/variants.vcf {output.vcf} 2>> {log}
        
        # Clean up
        rm -rf $TMPDIR
        """

rule run_pbsv:
    """
    Run PBSV structural variant caller.
    Specialized for PacBio long-read data.
    """
    input:
        bam = "results/mapped/{sample}.bam",
        bai = "results/mapped/{sample}.bam.bai",
        ref = config["reference"]["fasta"]
    output:
        vcf = "results/sv_calls/pbsv/{sample}.vcf",
        svsig = temp("results/sv_calls/pbsv/{sample}.svsig.gz")
    params:
        tr_arg = lambda wildcards: get_tr_argument(wildcards, "pbsv"),
        min_read_length = config["pbsv"]["min_read_length"],
        min_confidence = config["pbsv"]["min_confidence"],
        min_ref_span = config["pbsv"]["min_ref_span"],
        min_sv_length = config["pbsv"]["min_sv_length"]
    conda:
        "../../envs/environment.full.yaml"
    log:
        "logs/pbsv/{sample}.log"
    threads: config["threads"]
    shell:
        """
        # Discover signatures
        pbsv discover \
            {params.tr_arg} \
            --min-read-length {params.min_read_length} \
            {input.bam} \
            {output.svsig} \
            2> {log}
            
        # Call variants
        pbsv call \
            -j {threads} \
            --min-confidence {params.min_confidence} \
            --min-ref-span {params.min_ref_span} \
            --min-sv-length {params.min_sv_length} \
            {input.ref} \
            {output.svsig} \
            {output.vcf} \
            2>> {log}
        """

rule run_svdss:
    """
    Run SVDSS structural variant caller.
    Specialized for detecting structural variants using split read analysis.
    """
    input:
        bam = "results/mapped/{sample}.bam",
        bai = "results/mapped/{sample}.bam.bai",
        ref = config["reference"]["fasta"]
    output:
        vcf = "results/sv_calls/svdss/{sample}.vcf",
        smooth_bam = temp("results/sv_calls/svdss/{sample}.smooth.bam"),
        smooth_bai = temp("results/sv_calls/svdss/{sample}.smooth.bam.bai"),
        specifics = temp("results/sv_calls/svdss/{sample}.specifics.txt"),
        index = temp("results/sv_calls/svdss/{sample}.index")
    params:
        smooth_threads = config["svdss"]["smooth"]["threads"],
        search_threads = config["svdss"]["search"]["threads"],
        min_support = config["svdss"]["call"]["min_support"],
        min_size = config["svdss"]["call"]["min_size"]
    conda:
        "../../envs/environment.full.yaml"
    log:
        "logs/svdss/{sample}.log"
    threads: config["threads"]
    shell:
        """
        # Step 1: Smooth BAM file
        SVDSS smooth \
            --threads {params.smooth_threads} \
            --bam {input.bam} \
            --reference {input.ref} \
            > {output.smooth_bam} \
            2> {log}
            
        # Step 2: Index smoothed BAM
        samtools index {output.smooth_bam} {output.smooth_bai}
        
        # Step 3: Create index
        SVDSS index \
            --reference {input.ref} \
            --index {output.index} \
            2>> {log}
        
        # Step 4: Search for SVs
        SVDSS search \
            --threads {params.search_threads} \
            --index {output.index} \
            --bam {output.smooth_bam} \
            > {output.specifics} \
            2>> {log}
            
        # Step 5: Call variants
        SVDSS call \
            --reference {input.ref} \
            --bam {output.smooth_bam} \
            --sfs {output.specifics} \
            --min-cluster-weight {params.min_support} \
            --min-sv-length {params.min_size} \
            > {output.vcf} \
            2>> {log}
        """

rule run_debreak:
    """
    Run Debreak structural variant caller.
    Specialized for detecting structural variants in long-read data.
    """
    input:
        bam = "results/mapped/{sample}.bam",
        bai = "results/mapped/{sample}.bam.bai",
        ref = config["reference"]["fasta"]
    output:
        vcf = "results/sv_calls/debreak/{sample}.vcf"
    params:
        min_support = config["debreak"]["min_support"],
        min_size = config["debreak"]["min_size"],
        rescue_large_ins = config["debreak"]["rescue_large_ins"],
        rescue_dup = config["debreak"]["rescue_dup"],
        use_poa = config["debreak"]["use_poa"]
    container:
        config["containers"]["debreak"]
    log:
        "logs/debreak/{sample}.log"
    threads: config["threads"]
    shell:
        """
        # Run Debreak with configured parameters
        debreak \
            --min_support {params.min_support} \
            --min_size {params.min_size} \
            -t {threads} \
            --bam {input.bam} \
            -o ./ \
            {' --rescue_large_ins' if params.rescue_large_ins else ''} \
            {' --rescue_dup' if params.rescue_dup else ''} \
            {' --poa' if params.use_poa else ''} \
            --ref {input.ref} \
            2> {log}

        # Move VCF to final location
        mv debreak.vcf {output.vcf} 2>> {log}
            
        # Clean up temporary files
        rm -rf debreak-* \
              debreak_* \
              map_depth/ \
              sv_raw_calls/ \
              deletion-merged \
              duplication-merged \
              insertion-merged \
              inversion-merged \
              translocation-merged
        """
