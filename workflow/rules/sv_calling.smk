# SV calling rules for different callers

# Short-read callers
rule run_manta:
    input:
        bam = "results/mapped/{sample}.bam",
        ref = config["reference"]
    output:
        vcf = "results/sv_calls/manta/{sample}.vcf"
    shell:
        "configManta.py --bam {input.bam} --referenceFasta {input.ref} --runDir manta_tmp && \
         manta_tmp/runWorkflow.py"

rule run_delly:
    input:
        bam = "results/mapped/{sample}.bam",
        ref = config["reference"]
    output:
        vcf = "results/sv_calls/delly/{sample}.vcf"
    shell:
        "delly call -g {input.ref} {input.bam} > {output.vcf}"

# Add similar rules for lumpy and svaba

# Long-read callers
rule run_sniffles:
    input:
        bam = "results/mapped/{sample}.bam",
        ref = config["reference"]
    output:
        vcf = "results/sv_calls/sniffles/{sample}.vcf"
    shell:
        "sniffles --input {input.bam} --vcf {output.vcf}"

# Add similar rules for other long-read callers....... like SVDSS
