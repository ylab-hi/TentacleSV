# Rules for structural variant calling

# Helper function to get tandem repeats file path
def get_tandem_repeats():
    """Get tandem repeats file path if available."""
    return config["reference"].get("tandem_repeats", "")

# Helper function to construct tool-specific arguments
def get_tr_argument(wildcards, tool):
    """Get tool-specific tandem repeats argument."""
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
        config["containers"]["manta"]  # Use container instead of conda
    log:
        "logs/manta/{sample}.log"
    threads: config["threads"]
    shell:
        """
        # Create temporary directory for Manta
        TMPDIR=$(mktemp -d)
        
        # Configure Manta
        configManta.py \
            --bam {input.bam} \
            --referenceFasta {input.ref} \
            --runDir $TMPDIR \
            2> {log}
            
        # Run Manta workflow
        $TMPDIR/runWorkflow.py -j {threads} 2>> {log}
        
        # Move and clean up results
        mv $TMPDIR/results/variants/diploidSV.vcf.gz {output.vcf}.gz
        gunzip {output.vcf}.gz
        rm -rf $TMPDIR
        """

rule run_delly:
    input:
        bam = "results/mapped/{sample}.bam",
        bai = "results/mapped/{sample}.bam.bai",
        ref = config["reference"]["fasta"]
    output:
        vcf = "results/sv_calls/delly/{sample}.vcf"
    conda:
        "../../envs/environment.full.yaml"
    log:
        "logs/delly/{sample}.log"
    threads: config["threads"]
    shell:
        """
        # Run DELLY
        delly call \
            -g {input.ref} \
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
    container:
        config["containers"]["lumpy"]  # Use container instead of conda
    log:
        "logs/lumpy/{sample}.log"
    threads: config["threads"]
    shell:
        """
        lumpyexpress \
            -B {input.bam} \
            -o {output.vcf} \
            2> {log}
        """

rule run_svaba:
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
        # Run SVABA
        svaba run \
            -t {input.bam} \
            -p {threads} \
            -a svaba_out \
            -G {input.ref} \
            2> {log}
            
        # Move output to final location
        mv svaba_out.sv.vcf {output.vcf}
        """

# Long-read callers
rule run_cutesv:
    input:
        bam = "results/mapped/{sample}.bam",
        bai = "results/mapped/{sample}.bam.bai",
        ref = config["reference"]["fasta"]
    output:
        vcf = "results/sv_calls/cutesv/{sample}.vcf"
    conda:
        "../../envs/environment.full.yaml"
    log:
        "logs/cutesv/{sample}.log"
    threads: config["threads"]
    shell:
        """
        cuteSV \
            --threads {threads} \
            --genotype \
            --report_readid \
            -l 50 \
            -L 5000000 \
            -r 1000 \
            -q 20 \
            -s 3 \
            --max_cluster_bias_INS 1000 \
            --diff_ratio_merging_INS 0.9 \
            --max_cluster_bias_DEL 1000 \
            --diff_ratio_merging_DEL 0.5 \
            {input.bam} \
            {input.ref} \
            {output.vcf} \
            ./ \
            2> {log}
        """

rule run_sniffles:
    input:
        bam = "results/mapped/{sample}.bam",
        bai = "results/mapped/{sample}.bam.bai",
        ref = config["reference"]["fasta"]
    output:
        vcf = "results/sv_calls/sniffles/{sample}.vcf"
    params:
        tr_arg = lambda wildcards: get_tr_argument(wildcards, "sniffles")
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
            --minsupport auto \
            --minsvlen 50 \
            --threads {threads} \
            2> {log}
        """

rule run_svim:
    input:
        bam = "results/mapped/{sample}.bam",
        bai = "results/mapped/{sample}.bam.bai",
        ref = config["reference"]["fasta"]
    output:
        vcf = "results/sv_calls/svim/{sample}.vcf"
    conda:
        "../../envs/environment.full.yaml"
    log:
        "logs/svim/{sample}.log"
    shell:
        """
        # Create temporary directory for SVIM
        TMPDIR=$(mktemp -d)
        
        # Run SVIM alignment mode
        svim alignment \
            $TMPDIR \
            {input.bam} \
            {input.ref} \
            2> {log}
            
        # Move results
        mv $TMPDIR/variants.vcf {output.vcf} 2>> {log}
        
        # Clean up temporary directory
        rm -rf $TMPDIR
        """

rule run_pbsv:
    input:
        bam = "results/mapped/{sample}.bam",
        bai = "results/mapped/{sample}.bam.bai",
        ref = config["reference"]["fasta"]
    output:
        vcf = "results/sv_calls/pbsv/{sample}.vcf",
        svsig = temp("results/sv_calls/pbsv/{sample}.svsig.gz")
    params:
        tr_arg = lambda wildcards: get_tr_argument(wildcards, "pbsv")
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
            {input.bam} \
            {output.svsig} \
            2> {log}
            
        # Call variants
        pbsv call \
            -j {threads} \
            {input.ref} \
            {output.svsig} \
            {output.vcf} \
            2>> {log}
        """

rule run_svdss:
    input:
        bam = "results/mapped/{sample}.bam",
        bai = "results/mapped/{sample}.bam.bai",
        ref = config["reference"]["fasta"]
    output:
        vcf = "results/sv_calls/svdss/{sample}.vcf",
        smooth_bam = temp("results/sv_calls/svdss/{sample}.smooth.bam"),
        smooth_bai = temp("results/sv_calls/svdss/{sample}.smooth.bam.bai"),
        specifics = temp("results/sv_calls/svdss/{sample}.specifics.txt")
    conda:
        "../../envs/environment.full.yaml"
    log:
        "logs/svdss/{sample}.log"
    threads: config["threads"]
    shell:
        """
        # Step 1: Smooth BAM file
        SVDSS smooth \
            --threads {threads} \
            --bam {input.bam} \
            --reference {input.ref} > {output.smooth_bam} \
            2> {log}
            
        # Step 2: Index smoothed BAM
        samtools index {output.smooth_bam} {output.smooth_bai}
        
        # Step 3: Search for SVs
        SVDSS search \
            --threads {threads} \
            --index {input.ref} \
            --bam {output.smooth_bam} > {output.specifics} \
            2>> {log}
            
        # Step 4: Call variants
        SVDSS call \
            --reference {input.ref} \
            --bam {output.smooth_bam} \
            --sfs {output.specifics} > {output.vcf} \
            2>> {log}
        """

rule run_debreak:
    input:
        bam = "results/mapped/{sample}.bam",
        bai = "results/mapped/{sample}.bam.bai",
        ref = config["reference"]["fasta"]
    output:
        vcf = "results/sv_calls/debreak/{sample}.vcf"
    container:
        config["containers"]["debreak"]  # Use container instead of conda
    log:
        "logs/debreak/{sample}.log"
    threads: config["threads"]
    shell:
        """
        # Run Debreak with optimized parameters
        debreak \
            --min_support 2 \
            -t {threads} \
            --bam {input.bam} \
            -o ./ \
            --rescue_large_ins \
            --rescue_dup \
            --poa \
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
