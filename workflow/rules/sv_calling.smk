# Rules for structural variant calling

# Short-read callers
rule run_manta:
    input:
        bam = "results/mapped/{sample}.bam",
        bai = "results/mapped/{sample}.bam.bai",
        ref = config["reference"]["fasta"]
    output:
        vcf = "results/sv_calls/manta/{sample}.vcf"
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
        mv $TMPDIR/results/variants/tumorSV.vcf.gz {output.vcf}.gz
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
    log:
        "logs/lumpy/{sample}.log"
    threads: config["threads"]
    shell:
        """
        # Run LUMPY
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
            
        # Move and rename output
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
        ref = config["reference"]["fasta"],
        tr = config["reference"]["tandem_repeats"]
    output:
        vcf = "results/sv_calls/sniffles/{sample}.vcf"
    log:
        "logs/sniffles/{sample}.log"
    threads: config["threads"]
    shell:
        """
        sniffles \
            --input {input.bam} \
            --vcf {output.vcf} \
            --tandem-repeats {input.tr} \
            --reference {input.ref} \
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
    log:
        "logs/svim/{sample}.log"
    shell:
        """
        TMPDIR=$(mktemp -d)
        
        svim alignment \
            $TMPDIR \
            {input.bam} \
            {input.ref} \
            2> {log}
            
        # Filter and move results
        filter_vcf_based_on_quality.py \
            $TMPDIR/variants.vcf \
            10 > {output.vcf} \
            2>> {log}
            
        rm -rf $TMPDIR
        """

rule run_pbsv:
    input:
        bam = "results/mapped/{sample}.bam",
        bai = "results/mapped/{sample}.bam.bai",
        ref = config["reference"]["fasta"],
        tr = config["reference"]["tandem_repeats"]
    output:
        vcf = "results/sv_calls/pbsv/{sample}.vcf",
        svsig = temp("results/sv_calls/pbsv/{sample}.svsig.gz")
    log:
        "logs/pbsv/{sample}.log"
    threads: config["threads"]
    shell:
        """
        # Discover signatures
        pbsv discover \
            --tandem-repeats {input.tr} \
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
    log:
        "logs/svdss/{sample}.log"
    threads: config["threads"]
    shell:
        """
        # Smooth BAM
        SVDSS smooth \
            --threads {threads} \
            --bam {input.bam} \
            --reference {input.ref} > {output.smooth_bam} \
            2> {log}
            
        # Index smoothed BAM
        samtools index {output.smooth_bam} {output.smooth_bai}
        
        # Search for SVs
        SVDSS search \
            --threads {threads} \
            --index {input.ref}.fmd \
            --bam {output.smooth_bam} > {output.specifics} \
            2>> {log}
            
        # Call variants
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
    log:
        "logs/debreak/{sample}.log"
    threads: config["threads"]
    shell:
        """
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
