# TentacleSV

A Snakemake-based automated workflow for comprehensive structural variant (SV) analysis, built upon the powerful OctopusV toolkit. It streamlines the process of generating high-confidence structural variant call sets by integrating multiple SV callers and leveraging OctopusV's advanced BND correction and flexible merging capabilities.

TentacleSV provides:

* **Versatile Sequencing Support** - Handles both short-read and long-read sequencing data, supporting direct FASTQ processing or pre-aligned BAM files

* **Integrated Multi-caller Approach** - Incorporates major SV callers including Manta, Delly, Lumpy, SVaba for short reads, and CuteSV, PBSV, Sniffles2, SVIM, Debreak, SVDSS for long reads

* **Advanced SV Processing** - Features automated BND correction through OctopusV, standardization of SV annotations, and flexible merging strategies with customizable support thresholds

## Workflow Overview

TentacleSV implements a comprehensive pipeline that:

1. Performs read alignment using BWA-MEM2 (short reads) or Minimap2 (long reads)
2. Executes multiple SV callers in parallel
3. Corrects breakend (BND) annotations using OctopusV
4. Merges results with user-defined support criteria
5. Generates a high-confidence SV call set

All configurations are managed through a single YAML file, making it both user-friendly and highly customizable.
