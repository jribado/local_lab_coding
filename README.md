# Metagenomics QC pipeline

Authors: Jessica Ribado and Eli Moss
Aim: A simple wrapper for metagenomics QC using paired end reads. To use this pipeline, edit directories and program parameters in the config.json, and specify the proper path to config file in the submission script.

This program runs under the assumption samples are named (can be edited if other formats are more common):
PREFIX_R1.fastq.gz and PREFIX_R2.fastq.gz.

This script will create the following folders:
/PROJECT_DIR/qc/snakemake_files/logs
/PROJECT_DIR/qc/snakemake_files/benchmarks
/PROJECT_DIR/qc/00_qc_reports/pre_fastqc
/PROJECT_DIR/qc/00_qc_reports/post_fastqc
/PROJECT_DIR/qc/01_trimmed
/PROJECT_DIR/qc/02_dereplicate
/PROJECT_DIR/qc/03_interleave

Pipeline:
FASTQC -> TrimGalore -> SuperDeduper -> FASTQC/Interleave 

config.json parameters:
TrimGalore (Trim)
  - adaptors: Can specify either --illumina, --nextera or --small_rna as built in; alternative --adaptor [SEQUENCE] --adaptor2 [SEQUENCE]
  - quality: Minimun Phred score quality
  
SuperDeduper (Dereplicate)
  - start_trim: Starting location of or the starting base pair of the unique sequence ID. Program default 10.
  - unique_length: The length of the base pairs in the unique sequence ID (the number. of base pairs in the unique ID). Program default 25.
  
 FLASh (Interleave) 
  - min_length: Minimum required overlap length between two reads to provide a confident overlap. Default: 10bp.
  - max_length: Maximum overlap length expected in approximately 90% of read pairs. It is by default set to 70bp, which works well for 100bp reads generated from 180bp library (normal distribution of fragment lengths is assumed). Overlaps longer than maxOverlap are still considered as good overlaps, but the mismatch ratio (explained below) is calculated over the maxOverlap rather than the true overlap length. If you enter a value for maxOverlap, then the read length, fragment length and standard deviation of fragment lengths that you enter will be ignored for calculation of maxOverlap parameter. Default: 70bp.
  - mismatch: Maximum allowed ratio of the number of mismatches and the overlap length. An overlap with mismatch ratio higher than the set value is considered incorrect overlap and mates will not be merged. Any occurence of an "N" in any read is ignored and not counted towards the mismatches or overlap length. Default: 0.25.
