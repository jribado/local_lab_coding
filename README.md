# Metagenomics QC pipeline

Authors: Jessica Ribado and Eli Moss
Aim: A simple wrapper for metagenomics QC using paired end reads. To use this pipeline, edit directories and program parameters in the config.json, and specify the proper path to config file in the submission script.

This program runs under the assumption samples are named (can be edited if other formats are more common):
PREFIX_R1.fastq.gz and PREFIX_R2.fastq.gz.

This script will create the following folders:
/PROJECT_DIR/qc/00_qc_reports/pre_fastqc
/PROJECT_DIR/qc/00_qc_reports/post_fastqc
/PROJECT_DIR/qc/01_trimmed
/PROJECT_DIR/qc/02_dereplicate
/PROJECT_DIR/qc/03_interleave_seqtk

Pipeline:
FASTQC -> Trim low quality bases -> PCR Dereplication -> FASTQC/Interleave

Please note that SuperDeduper could not be wrapped into the

config.json parameters:
Working directory (Output files path)
Raw fastq directory (Input file path - will parse names from directory)

TrimGalore (Trim low quality bases)
  - adaptors: Can specify either --illumina, --nextera or --small_rna as built in; alternative --adaptor [SEQUENCE] --adaptor2 [SEQUENCE]
  - quality: Minimun Phred score quality

Seqkit (Dereplicate PCR amplicons)
  - start_trim: Starting location of or the starting base pair of the unique sequence ID.

Seqtk (Interleave)
  - no parameters
