import os,re

'''
Authors: Jessica Ribado and Eli Moss
Aim: A simple wrapper for metagenomics QC using paired end reads. To use this pipeline, edit parameters in the config.yaml, and specify the proper path to config file in the submission script.

This program runs under the assumption samples are named:
PREFIX_R1.fastq.gz and PREFIX_R2.fastq.gz.

This script will create the following folders:
PROJECT_DIR/qc/snakemake_files/logs
PROJECT_DIR/qc/snakemake_files/benchmarks
PROJECT_DIR/qc/00_qc_reports/pre_fastqc
PROJECT_DIR/qc/00_qc_reports/post_fastqc
PROJECT_DIR/qc/01_trimmed
PROJECT_DIR/qc/02_dereplicate
PROJECT_DIR/qc/03_interleave

Run:
snakemake --jobs 100 \
		  -latency-wait 60 \
		  --configfile CONFIG_PATH/config.json \
		  --cluster-config CONFIG_PATH/clusterconfig.json \
		  --profile scg	
'''

################################################################################
# specify project directories
DATA_DIR	= config["raw_reads_directory"]
PROJECT_DIR = config["output_directory"]
LOGS_DIR	= os.path.join(PROJECT_DIR, "qc/snakemake_files/logs")
BENCH_DIR   = os.path.join(PROJECT_DIR, "qc/snakemake_files/benchmarks")
DEDUP_DIR   = "/srv/gsfs0/projects/bhatt/tools/moss_tools/qc/Super-Deduper"


################################################################################
# get the names of the files in the directory
FILES = [f for f in os.listdir(DATA_DIR) if f.endswith('fastq.gz')]
SAMPLE_PREFIX = list(set([i.split('_R', 1)[0] for i in FILES]))

################################################################################
# specify which rules do not need to be submitted to the cluster
localrules: multiQC_report, interleave

rule all:
	input:
		expand(os.path.join(PROJECT_DIR, "qc/00_qc_reports/{qc}_multiqc_report.html"), qc=['pre', 'post']),
		expand(os.path.join(PROJECT_DIR + "/qc/03_interleave/{sample}_trimmed_{mismatch}.extendedFrags.fastq.gz"), sample=SAMPLE_PREFIX, mismatch=config['interleave']['mismatch'])

################################################################################
rule pre_fastqc:
	input: os.path.join(DATA_DIR, "{sample}_R{read}.fastq.gz")
	output: os.path.join(PROJECT_DIR,  "qc/00_qc_reports/pre_fastqc/{sample}_R{read}_fastqc.html")
	threads: 1
	log: os.path.join(LOGS_DIR, "pre_fastqc_{sample}_R{read}")
	shell: """
	   mkdir -p {PROJECT_DIR}/qc/00_qc_reports/pre_fastqc/
	   module load java/latest
	   module load fastqc/0.11.2
	   fastqc {input} --outdir {PROJECT_DIR}/qc/00_qc_reports/pre_fastqc/
	"""

################################################################################
rule trim_galore:
	input:
		fwd = os.path.join(DATA_DIR, "{sample}_R1.fastq.gz"),
		rev = os.path.join(DATA_DIR, "{sample}_R2.fastq.gz")
	output:
		fwd = os.path.join(PROJECT_DIR, "qc/01_trimmed/{sample}_R1_val_1.fq.gz"),
		rev = os.path.join(PROJECT_DIR, "qc/01_trimmed/{sample}_R2_val_2.fq.gz")
	threads: 4
	params:
		adaptor = config['trim_galore']['adaptors'],
		q_min   = config['trim_galore']['quality']
	log: os.path.join(LOGS_DIR, "trimGalore_{sample}")
	#benchmark: os.path.join(BENCH_DIR, "trimGalore_{sample}.txt")
	shell: """
		 mkdir -p {PROJECT_DIR}/qc/01_trimmed/
		 module load fastqc/0.11.2
		 module load trim_galore/0.4.2
		 trim_galore --{params.adaptor} \
					 --quality {params.q_min} \
					 --output_dir {PROJECT_DIR}/qc/01_trimmed/ \
					 --paired {input.fwd} {input.rev}
		 """

################################################################################
rule dereplicate:
	input:
	 	fwd = rules.trim_galore.output.fwd,
		rev = rules.trim_galore.output.rev
	output:
	 	fwd = os.path.join(PROJECT_DIR, "qc/02_dereplicate/{sample}_nodup_PE1.fastq"),
		rev = os.path.join(PROJECT_DIR, "qc/02_dereplicate/{sample}_nodup_PE2.fastq")
	threads: 2
	params:
		start	= config['dereplicate']['start_trim'],
		length	= config['dereplicate']['unique_length'],
		prefix	= "{sample}"
	log: os.path.join(LOGS_DIR, "derep_{sample}")
	#benchmark: os.path.join(BENCH_DIR, "superDeduper_{sample}.txt")
	shell: """
		{DEDUP_DIR}/super_deduper -1 {input.fwd} -2 {input.rev} \
			-p {PROJECT_DIR}/qc/02_dereplicate/{params.prefix} \
			-- start {params.start} \
			-- length {params.length}
		 """

################################################################################
rule post_fastqc:
	input:  rules.dereplicate.output
	output: os.path.join(PROJECT_DIR,  "qc/00_qc_reports/post_fastqc/{sample}_nodup_PE{read}_fastqc.html")
	threads: 4
	log: os.path.join(LOGS_DIR + "post_fastqc_{sample}_R{read}")
	shell: """
	   mkdir -p {PROJECT_DIR}/qc/00_qc_reports/post_fastqc/
	   module load fastqc/0.11.2
	   fastqc {input} -f fastq --outdir {PROJECT_DIR}/qc/00_qc_reports/post_fastqc/
	  """

################################################################################
rule multiQC_report:
	input:  os.path.join(PROJECT_DIR, "qc/00_qc_reports/{qc}_fastqc/")
	output: os.path.join(PROJECT_DIR, "qc/00_qc_reports/{qc}_multiqc_report.html")
	threads: 1
	shell:
		"""
		module load multiqc/1.5
		multiqc -f {input} -o {PROJECT_DIR}/qc/00_qc_reports/ -n {wildcards.qc}_multiqc_report.html
		"""

################################################################################
rule interleave:
	input:  rules.dereplicate.output
	output: os.path.join(PROJECT_DIR,  "qc/03_interleave/{sample}_trimmed_{mismatch}.extendedFrags.fastq.gz")
	params:
		min_len   = config['interleave']['min_length'],
		max_len   = config['interleave']['max_length'],
		mismatch  = config['interleave']['mismatch'],
		prefix	  = "{sample}"
	log: os.path.join(LOGS_DIR, "flash_{mismatch}_mergepairs.log")
	#benchmark: os.path.join(BENCH_DIR, "flash_{sample}.txt")
	shell: """
		mkdir -p {PROJECT_DIR}/qc/03_interleave/
		module load flash/1.2.11
		flash -m {params.min_len} -M {params.max_len} \
			  -x {params.mismatch} \
			  -o {params.prefix}_trimmed_{params.mismatch} \
			  -d {PROJECT_DIR}/qc/03_interleave -z {input}
	"""
