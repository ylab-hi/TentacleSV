# TentacleSV

![PyPI](https://img.shields.io/badge/pypi-v1.0.0-blue)
![Snakemake](https://img.shields.io/badge/snakemake-≥6.0.0-brightgreen.svg)
![License](https://img.shields.io/badge/license-Apache--2.0-green)
![Python](https://img.shields.io/badge/python-3.10%20|%203.11%20|%203.12-blue)
![Platform](https://img.shields.io/badge/platform-linux%20|%20osx-lightgrey)
![GitHub Workflow Status](https://img.shields.io/badge/CI-passing-brightgreen)

A Snakemake-based automated workflow for generating high-confidence structural variant (SV) call sets, powered by [OctopusV](https://github.com/ylab-hi/octopusV). TentacleSV streamlines multi-caller SV analysis by integrating state-of-the-art callers with OctopusV's advanced BND correction and flexible merging capabilities.

## 🌟 Key Features

- **Universal Input Support**: Process both FASTQ files and pre-aligned BAM files
- **Multi-platform Analysis**: Comprehensive support for short-read and long-read sequencing data
- **Automated Pipeline**: From raw sequencing data to high-confidence SV call sets
- **Powered by OctopusV**: Advanced BND correction and flexible merging strategies

## 🔧 Supported Tools

**Short-read Callers**:
- Manta
- Delly
- Lumpy
- SVaba

**Long-read Callers**:
- CuteSV
- PBSV
- Sniffles2
- SVIM
- Debreak
- SVDSS

## 📦 Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/TentacleSV.git
cd TentacleSV

# Create conda environment
mamba env create -f environment.yaml
conda activate tentaclesv
```

## 🚀 Quick Start

1. **Configure Your Analysis**

Edit `config/config.yaml` to specify your input data and parameters:

```yaml
# For short-read data
seq_type: "short"
samples:
  sample1:
    type: "fastq"
    fq1: "path/to/R1.fastq.gz"
    fq2: "path/to/R2.fastq.gz"

# For long-read data
seq_type: "long_pacbio"  # or "long_ont"
samples:
  sample1:
    type: "fastq"
    fq1: "path/to/reads.fastq.gz"
```

2. **Run the Pipeline**

```bash
snakemake --cores 8 \
          --use-conda \
          --use-singularity \
          --singularity-args "--bind /projects,/home" \
          --conda-frontend mamba \
          --rerun-incomplete \
          --keep-going \
          -s workflow/Snakefile
```

## 📊 Output

TentacleSV generates:
- High-confidence SV call set (`results/merged/{sample}.vcf`)
- Quality metrics and logs
- Optional UpSet plots showing caller overlap

## 📝 Citation

If you use TentacleSV in your research, please cite:

```bibtex
@software{TentacleSV2024,
  author = {Qing Gao},
  title = {TentacleSV: A Comprehensive SV Analysis Pipeline},
  year = {2024},
  publisher = {GitHub},
  url = {https://github.com/ylab-hi/TentacleSV}
}
```

## 📄 License

This project is licensed under the Apache 2.0 License - see the [LICENSE](LICENSE) file for details.

## 🤝 Contributing

Contributions are welcome! Please feel free to submit a Pull Request.
