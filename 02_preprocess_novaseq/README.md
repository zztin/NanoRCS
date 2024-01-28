# Preprocessing: 02_preprocessing_novaseq

## Introduction
Illumina NovaSeq paired-end sequencing files are processed to produce mapped bam files, read length distribution, and other statistics in this folder with a Snankemake pipeline. Example output is directed to `NanoRCS/output/processed_data/novaseq`. Downstream analysis are referred to the relative path from this folder.

An example submission command: `snakemake --snakefile Snakemake --configfiles config-novaseq-hs37d5.yaml --profile slurm`

Here, `--profile slurm` is used. The related snakemake profile could be found at top level folder: `NanoRCS/snakemake-profile`. More snakemake related documentations please be referred to [Snakemake Documentation](https://snakemake.readthedocs.io/en/stable/).

