# cfDNA Analyses: 03_snv_error_rate

## Introduction
Analysis is related to: Method paragraph "cfDNA SNV detection error rate" in the paper.

Cell-free DNA of three HCs with known genome references were analysed. All reads overlapping SNP in any healthy controls is excluded. 

SNV error rate is derived from counting CIGAR differences of all remaining reads on chromosome 1-22 (in aggregated batch of 100,000 reads, or per read).

SNV error rate for 2 different sequencing techniques, with 4 different way of data processing are dervied. 

Output from analysis `01_preprocess_nanorcs`,`02_preprocess_novaseq` are required. 

# Steps
* Step 1: Generate BAM and VCF files from HCs. Filter VCF files.
* Step 2: Overlap BAM files with known VCF file.
  * Output: results/overlap/{sample}_overlap_all_exclude.bam
* Step 3: Convert BAM files to pickle files that contains curated information per read.
* Step 4: get_error_rate.py: generate source data
* Step 5: plotting (in folder 10_figures)

##  Step 1: Generate BAM and VCF files from HCs. Filter VCF files.

1. Sarek3 is used to generate VCF files. 
2. Filtering variants and convert to BED files: `XXXX`. 

## Step 2: Get first alignment in NanoRCS to resemble raw nanopore sequencing (Named Raw NanoRCS)
1. Use first_align_nanorcs_wrapper.sh to submit first_align_nanorcs.py per chromosome.
2. Slurm parallel submission option: check first_align_nanorcs_wrapper.sh and first_align_nanorcs_parallel_submission_exmaple.sh 

## Step 3: Overlap BAM files with known VCF file.
1. Modify config file to include relevant BAM file paths and BED file path.
2. `snakemake --Snakefile Snakemake --cores all --configfiles configs/<config-template.yaml>` More snakemake related documentations please be referred to [Snakemake Documentation](https://snakemake.readthedocs.io/en/stable/).

## Step 4: Convert BAM files to pickle files that contains curated information per read.
`error_rate_bam_to_pickle.py`

Execute environment: 
##  Step 5: get_error_rate.py: generate source data
`error_rate_pickle_to_csv.py`
Execute environment: 
