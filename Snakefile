import os
import re
import sys

configfile: "config.yaml"
gatk = config["gatk_path"]
ref_path0 = config["reference_panel_path"]
def fa_dict(path0: str) -> str:
	dir_path = os.path.dirname(path0)
	file_name = os.path.splitext(os.path.basename(path0))[0]
	output_file = os.path.join(dir_path, file_name + ".dict")
	return output_file

ref_dict_path = fa_dict(ref_path0)

Trio_ID = {config['CHILD_R1_path']}.split("/")[-1].split("_")[0]

db = ""
for i in config["category"].split("-"):
    db += f"-{i} Y "

if ('ASD' in {config['panel_type']}):
	getTrioVar_flag = 1
else:
	getTrioVar_flag = 0

rule end:
	input:

# ====================== fastp ======================
rule FUO_fastp:
	input:
		f"{config['FUO_R1_path']}",
		f"{config['FUO_R2_path']}"
	output:
		f"{config['output_dir']}/fastp/{config['FUO_sample_name']}.R1.fastp.gz",
		f"{config['output_dir']}/fastp/{config['FUO_sample_name']}.R2.fastp.gz",
		f"{config['output_dir']}/fastp/{config['FUO_sample_name']}_fastp.html",
		f"{config['output_dir']}/fastp/{config['FUO_sample_name']}_fastp.json"
	conda:
		"./first_step_mamba.yml"
	shell:
		f"""fastp -i {config['FUO_R1_path']} -I {config['FUO_R2_path']} \\
		-o {config['output_dir']}/fastp/{config['FUO_sample_name']}.R1.fastp.gz -O {config['output_dir']}/fastp/{config['FUO_sample_name']}.R2.fastp.gz \\
		-h {config['output_dir']}/fastp/{config['FUO_sample_name']}_fastp.html -j {config['output_dir']}/fastp/{config['FUO_sample_name']}_fastp.json"""

rule MU0_fastp:
	input:
		f"{config['MU0_R1_path']}",
		f"{config['MU0_R2_path']}"
	output:
		f"{config['output_dir']}/fastp/{config['MU0_sample_name']}.R1.fastp.gz",
		f"{config['output_dir']}/fastp/{config['MU0_sample_name']}.R2.fastp.gz",
		f"{config['output_dir']}/fastp/{config['MU0_sample_name']}_fastp.html",
		f"{config['output_dir']}/fastp/{config['MU0_sample_name']}_fastp.json"
	conda:
		"./first_step_mamba.yml"
	shell:
		f"""fastp -i {config['MU0_R1_path']} -I {config['MU0_R2_path']} \\
		-o {config['output_dir']}/fastp/{config['MU0_sample_name']}.R1.fastp.gz -O {config['output_dir']}/fastp/{config['MU0_sample_name']}.R2.fastp.gz \\
		-h {config['output_dir']}/fastp/{config['MU0_sample_name']}_fastp.html -j {config['output_dir']}/fastp/{config['MU0_sample_name']}_fastp.json"""

rule Child_fastp:
	input:
		f"{config['CHILD_R1_path']}",
		f"{config['CHILD_R2_path']}"
	output:
		f"{config['output_dir']}/fastp/{config['CHILD_sample_name']}.R1.fastp.gz",
		f"{config['output_dir']}/fastp/{config['CHILD_sample_name']}.R2.fastp.gz",
		f"{config['output_dir']}/fastp/{config['CHILD_sample_name']}_fastp.html",
		f"{config['output_dir']}/fastp/{config['CHILD_sample_name']}_fastp.json"
	conda:
		"./first_step_mamba.yml"
	shell:
		f"""fastp -i {config['CHILD_R1_path']} -I {config['CHILD_R2_path']} \\
		-o {config['output_dir']}/fastp/{config['CHILD_sample_name']}.R1.fastp.gz -O {config['output_dir']}/fastp/{config['CHILD_sample_name']}.R2.fastp.gz \\
		-h {config['output_dir']}/fastp/{config['CHILD_sample_name']}_fastp.html -j {config['output_dir']}/fastp/{config['CHILD_sample_name']}_fastp.json"""

# ====================== fastQC ======================
rule FUO_fastQC:
	input:
		f"{config['output_dir']}/fastp/{config['FUO_sample_name']}.R1.fastp.gz",
		f"{config['output_dir']}/fastp/{config['FUO_sample_name']}.R2.fastp.gz"
	output:
		f"{config['output_dir']}/fastqc/{config['FUO_sample_name']}_fastqc.html",
		f"{config['output_dir']}/fastqc/{config['FUO_sample_name']}_fastqc.zip"
	conda:
		"./first_step_mamba.yml"
	shell:
		f"""fastqc {config['output_dir']}/fastp/{config['FUO_sample_name']}.R1.fastp.gz {config['output_dir']}/fastp/{config['FUO_sample_name']}.R2.fastp.gz \\
		-o {config['output_dir']}/fastqc/"""

rule MU0_fastQC:
	input:
		f"{config['output_dir']}/fastp/{config['MU0_sample_name']}.R1.fastp.gz",
		f"{config['output_dir']}/fastp/{config['MU0_sample_name']}.R2.fastp.gz"
	output:
		f"{config['output_dir']}/fastqc/{config['MU0_sample_name']}_fastqc.html",
		f"{config['output_dir']}/fastqc/{config['MU0_sample_name']}_fastqc.zip"
	conda:
		"./first_step_mamba.yml"
	shell:
		f"""fastqc {config['output_dir']}/fastp/{config['MU0_sample_name']}.R1.fastp.gz {config['output_dir']}/fastp/{config['MU0_sample_name']}.R2.fastp.gz \\
		-o {config['output_dir']}/fastqc/"""

rule Child_fastQC:
	input:
		f"{config['output_dir']}/fastp/{config['CHILD_sample_name']}.R1.fastp.gz",
		f"{config['output_dir']}/fastp/{config['CHILD_sample_name']}.R2.fastp.gz"
	output:
		f"{config['output_dir']}/fastqc/{config['CHILD_sample_name']}_fastqc.html",
		f"{config['output_dir']}/fastqc/{config['CHILD_sample_name']}_fastqc.zip"
	conda:
		"./first_step_mamba.yml"
	shell:
		f"""fastqc {config['output_dir']}/fastp/{config['CHILD_sample_name']}.R1.fastp.gz {config['output_dir']}/fastp/{config['CHILD_sample_name']}.R2.fastp.gz \\
		-o {config['output_dir']}/fastqc/"""

# ====================== BWA ======================
rule FUO_align_readers:
	input:
		f"{config['reference_panel_path']}.bwt",
		f"{config['reference_panel_path']}.pac",
		f"{config['reference_panel_path']}.ann",
		f"{config['reference_panel_path']}.sa",
		f"{config['reference_panel_path']}.amb",
		f"{config['reference_panel_path']}",
		f"{config['output_dir']}/fastp/{config['FUO_sample_name']}.R1.fastp.gz",
		f"{config['output_dir']}/fastp/{config['FUO_sample_name']}.R2.fastp.gz"
	output:
		f"{config['output_dir']}/alignments/{config['FUO_sample_name']}.bwa.bam"
	conda:
		"./first_step_mamba.yml"
	threads: 16
	shell:
		f"bwa mem -t {threads}" + f" {config['reference_panel_path']} {config['output_dir']}/fastp/{config['FUO_sample_name']}.R1.fastp.gz {config['output_dir']}/fastp/{config['FUO_sample_name']}.R2.fastp.gz | samtools sort > {config['output_dir']}/alignments/{config['FUO_sample_name']}.bwa.bam" # bwa mem will align the reads to the reference panel

rule MU0_align_readers:
	input:
		f"{config['reference_panel_path']}.bwt",
		f"{config['reference_panel_path']}.pac",
		f"{config['reference_panel_path']}.ann",
		f"{config['reference_panel_path']}.sa",
		f"{config['reference_panel_path']}.amb",
		f"{config['reference_panel_path']}",
		f"{config['output_dir']}/fastp/{config['MU0_sample_name']}.R1.fastp.gz",
		f"{config['output_dir']}/fastp/{config['MU0_sample_name']}.R2.fastp.gz"
	output:
		f"{config['output_dir']}/alignments/{config['MU0_sample_name']}.bwa.bam"
	conda:
		"./first_step_mamba.yml"
	threads: 16
	shell:
		f"bwa mem -t {threads}" + f" {config['reference_panel_path']} {config['output_dir']}/fastp/{config['MU0_sample_name']}.R1.fastp.gz {config['output_dir']}/fastp/{config['MU0_sample_name']}.R2.fastp.gz | samtools sort > {config['output_dir']}/alignments/{config['MU0_sample_name']}.bwa.bam" # bwa mem will align the reads to the reference panel

rule Child_align_readers:
	input:
		f"{config['reference_panel_path']}.bwt",
		f"{config['reference_panel_path']}.pac",
		f"{config['reference_panel_path']}.ann",
		f"{config['reference_panel_path']}.sa",
		f"{config['reference_panel_path']}.amb",
		f"{config['reference_panel_path']}",
		f"{config['output_dir']}/fastp/{config['CHILD_sample_name']}.R1.fastp.gz",
		f"{config['output_dir']}/fastp/{config['CHILD_sample_name']}.R2.fastp.gz"
	output:
		f"{config['output_dir']}/alignments/{config['CHILD_sample_name']}.bwa.bam"
	conda:
		"./first_step_mamba.yml"
	threads: 16
	shell:
		f"bwa mem -t {threads}" + f" {config['reference_panel_path']} {config['output_dir']}/fastp/{config['CHILD_sample_name']}.R1.fastp.gz {config['output_dir']}/fastp/{config['CHILD_sample_name']}.R2.fastp.gz | samtools sort > {config['output_dir']}/alignments/{config['CHILD_sample_name']}.bwa.bam" # bwa mem will align the reads to the reference panel
# ====================== AddOrReplaceReadGroups ======================
rule FUO_AddOrReplaceReadGroups: # Provide information for BQSR
	input:
		f"{config['output_dir']}/alignments/{config['FUO_sample_name']}.bwa.bam"
	output:
		f"{config['output_dir']}/alignments/{config['FUO_sample_name']}.bwa.rg.bam"
	shell:
		f"""{gatk} AddOrReplaceReadGroups \\
			I={config['output_dir']}/alignments/{config['FUO_sample_name']}.bwa.bam \\
			O={config['output_dir']}/alignments/{config['FUO_sample_name']}.bwa.rg.bam \\
			RGID={config['FUO_sample_name']} \\
			RGLB=lib1 \\
			RGPL={config['platform']} \\
			RGPU={config['lane']} \\
			RGSM={config['FUO_sample_name']} \\
			CREATE_INDEX=True
"""

rule MU0_AddOrReplaceReadGroups: # Provide information for BQSR
	input:
		f"{config['output_dir']}/alignments/{config['MU0_sample_name']}.bwa.bam"
	output:
		f"{config['output_dir']}/alignments/{config['MU0_sample_name']}.bwa.rg.bam"
	shell:
		f"""{gatk} AddOrReplaceReadGroups \\
			I={config['output_dir']}/alignments/{config['MU0_sample_name']}.bwa.bam \\
			O={config['output_dir']}/alignments/{config['MU0_sample_name']}.bwa.rg.bam \\
			RGID={config['MU0_sample_name']} \\
			RGLB=lib1 \\
			RGPL={config['platform']} \\
			RGPU={config['lane']} \\
			RGSM={config['MU0_sample_name']} \\
			CREATE_INDEX=True
"""

rule Child_AddOrReplaceReadGroups:
	input:
		f"{config['output_dir']}/alignments/{config['CHILD_sample_name']}.bwa.bam"
	output:
		f"{config['output_dir']}/alignments/{config['CHILD_sample_name']}.bwa.rg.bam"
	shell:
		f"""{gatk} AddOrReplaceReadGroups \\
			I={config['output_dir']}/alignments/{config['CHILD_sample_name']}.bwa.bam \\
			O={config['output_dir']}/alignments/{config['CHILD_sample_name']}.bwa.rg.bam \\
			RGID={config['CHILD_sample_name']} \\
			RGLB=lib1 \\
			RGPL={config['platform']} \\
			RGPU={config['lane']} \\
			RGSM={config['CHILD_sample_name']} \\
			CREATE_INDEX=True
"""

# ====================== MarkDuplicatesSpark ======================
rule FUO_mark_duplicates: # Notice MarkDuplicatesSpark is not BETA
	input:
		f"{config['output_dir']}/alignments/{config['FUO_sample_name']}.bwa.rg.bam"
	output:
		f"{config['output_dir']}/alignments/{config['FUO_sample_name']}.bwa.markdup.rg.bam",
		f"{config['output_dir']}/alignments/{config['FUO_sample_name']}.bwa.markdup.rg.metrics"
	shell:
		f"""{gatk} MarkDuplicatesSpark \\ 
			-I {config['output_dir']}/alignments/{config['FUO_sample_name']}.bwa.rg.bam \\
			-O {config['output_dir']}/alignments/{config['FUO_sample_name']}.bwa.markdup.rg.bam \\
			-M {config['output_dir']}/alignments/{config['FUO_sample_name']}.bwa.markdup.rg.metrics \\
			--create-output-bam-index true
""" # sambamba markdup will mark the duplicates in the bam file

			#--RGPU -PU: Read-Group platform unit (eg. run barcode)
			#--RGLB -LB: Read-Group library
			#--RGPL -PL: Read-Group platform (eg. illumina, solid)
			#--RGSM -SM: Read-Group sample name
			#--RGID -ID: Read-Group ID
			#--RGCN -CN: Read-Group sequencing center name
			#--RGDS -DS: Read-Group description
			#--RGDT -DT: Read-Group run date
			#--RGPI -PI: Read-Group predicted insert size
			#--RGPG -PG: Read-Group program group
			#--RGPM -PM: Read-Group platform model
			#--RGKS -KS: Read-Group key sequence
			#--RGFO -FO: Read-Group flow order

rule MU0_MarkDuplicatesSpark: # Notice MarkDuplicatesSpark is not BETA
	input:
		f"{config['output_dir']}/alignments/{config['MU0_sample_name']}.bwa.rg.bam"
	output:
		f"{config['output_dir']}/alignments/{config['MU0_sample_name']}.bwa.markdup.rg.bam",
		f"{config['output_dir']}/alignments/{config['MU0_sample_name']}.bwa.markdup.rg.metrics"
	shell:
		f"""{gatk} MarkDuplicatesSpark \\ 
			-I {config['output_dir']}/alignments/{config['MU0_sample_name']}.bwa.rg.bam \\
			-O {config['output_dir']}/alignments/{config['MU0_sample_name']}.bwa.markdup.rg.bam \\
			-M {config['output_dir']}/alignments/{config['MU0_sample_name']}.bwa.markdup.rg.metrics \\
			--create-output-bam-index true
""" 

rule Child_MarkDuplicatesSpark:
	input:
		f"{config['output_dir']}/alignments/{config['CHILD_sample_name']}.bwa.rg.bam"
	output:
		f"{config['output_dir']}/alignments/{config['CHILD_sample_name']}.bwa.markdup.rg.bam",
		f"{config['output_dir']}/alignments/{config['CHILD_sample_name']}.bwa.markdup.rg.metrics"
	shell:
		f"""{gatk} MarkDuplicatesSpark \\ 
			-I {config['output_dir']}/alignments/{config['CHILD_sample_name']}.bwa.rg.bam \\
			-O {config['output_dir']}/alignments/{config['CHILD_sample_name']}.bwa.markdup.rg.bam \\
			-M {config['output_dir']}/alignments/{config['CHILD_sample_name']}.bwa.markdup.rg.metrics \\
			--create-output-bam-index true
"""


# ====================== Index MarkDuplicates' BAM ======================
rule FUO_mark_duplicates_index_bam:
	input:
		f"{config['output_dir']}/alignments/{config['FUO_sample_name']}.bwa.markdup.rg.bam"
	output:
		f"{config['output_dir']}/alignments/{config['FUO_sample_name']}.bwa.markdup.rg.bam.bai"
	conda:
		"./first_step_mamba.yml"
	shell:
		f"samtools index {config['output_dir']}/alignments/{config['FUO_sample_name']}.bwa.markdup.rg.bam" 

rule MU0_MarkDuplicates_index_bam:
	input:
		f"{config['output_dir']}/alignments/{config['MU0_sample_name']}.bwa.markdup.rg.bam"
	output:
		f"{config['output_dir']}/alignments/{config['MU0_sample_name']}.bwa.markdup.rg.bam.bai"
	conda:
		"./first_step_mamba.yml"
	shell:
		f"samtools index {config['output_dir']}/alignments/{config['MU0_sample_name']}.bwa.markdup.rg.bam"

rule Child_MarkDuplicates_index_bam:
	input:
		f"{config['output_dir']}/alignments/{config['CHILD_sample_name']}.bwa.markdup.rg.bam"
	output:
		f"{config['output_dir']}/alignments/{config['CHILD_sample_name']}.bwa.markdup.rg.bam.bai"
	conda:
		"./first_step_mamba.yml"
	shell:
		f"samtools index {config['output_dir']}/alignments/{config['CHILD_sample_name']}.bwa.markdup.rg.bam"

## =========================== BQSR ============================
# ====================== BaseRecalibrator ======================
rule FUO_BaseRecalibrator:
	input:
		f"{config['output_dir']}/alignments/{config['FUO_sample_name']}.bwa.markdup.rg.bam",
		f"{config['output_dir']}/alignments/{config['FUO_sample_name']}.bwa.markdup.rg.bam.bai",
		f"{config['reference_panel_path']}",
		f"{config['reference_panel_path']}.fai",
		f"{config['known_sites']}",
		f"{config['known_sites']}.tbi"
	output:
		f"{config['output_dir']}/alignments/{config['FUO_sample_name']}.bwa.markdup.rg.bam.table"
	shell:
		f"""{gatk} BaseRecalibrator \\
			-R {config['reference_panel_path']} \\
			-I {config['output_dir']}/alignments/{config['FUO_sample_name']}.bwa.markdup.rg.bam \\
			-O {config['output_dir']}/alignments/{config['FUO_sample_name']}.bwa.markdup.rg.bam.table \\
			--use-original-qualities \\
			--known-sites {config['known_sites']} \\
			--known-sites {config['known_sites2']}""" # known-sites2 -> indel

rule MUO_BaseRecalibrator:
	input:
		f"{config['output_dir']}/alignments/{config['MU0_sample_name']}.bwa.markdup.rg.bam",
		f"{config['output_dir']}/alignments/{config['MU0_sample_name']}.bwa.markdup.rg.bam.bai",
		f"{config['reference_panel_path']}",
		f"{config['reference_panel_path']}.fai",
		f"{config['known_sites']}",
		f"{config['known_sites']}.tbi"
	output:
		f"{config['output_dir']}/alignments/{config['MU0_sample_name']}.bwa.markdup.rg.bam.table"
	shell:
		f"""{gatk} BaseRecalibrator \\
			-R {config['reference_panel_path']} \\
			-I {config['output_dir']}/alignments/{config['MU0_sample_name']}.bwa.markdup.rg.bam \\
			-O {config['output_dir']}/alignments/{config['MU0_sample_name']}.bwa.markdup.rg.bam.table \\
			--use-original-qualities \\
			--known-sites {config['known_sites']}"""

rule Child_BaseRecalibrator:
	input:	
		f"{config['output_dir']}/alignments/{config['CHILD_sample_name']}.bwa.markdup.rg.bam",
		f"{config['output_dir']}/alignments/{config['CHILD_sample_name']}.bwa.markdup.rg.bam.bai",
		f"{config['reference_panel_path']}",
		f"{config['reference_panel_path']}.fai",
		f"{config['known_sites']}",
		f"{config['known_sites']}.tbi"
	output:
		f"{config['output_dir']}/alignments/{config['CHILD_sample_name']}.bwa.markdup.rg.bam.table"
	shell:
		f"""{gatk} BaseRecalibrator \\
			-R {config['reference_panel_path']} \\
			-I {config['output_dir']}/alignments/{config['CHILD_sample_name']}.bwa.markdup.rg.bam \\
			-O {config['output_dir']}/alignments/{config['CHILD_sample_name']}.bwa.markdup.rg.bam.table \\
			--use-original-qualities \\
			--known-sites {config['known_sites']}"""
# ====================== ApplyBQSR ======================

rule FUO_ApplyBQSR:
	input:
		f"{config['output_dir']}/alignments/{config['FUO_sample_name']}.bwa.markdup.rg.bam.table"
	output:
		f"{config['output_dir']}/alignments/{config['FUO_sample_name']}.bwa.markdup.rg.bqsr.bam"
	shell:
		f"""{gatk} ApplyBQSRSpark -R {config['reference_panel_path']} \\
			-I {config['output_dir']}/alignments/{config['FUO_sample_name']}.bwa.markdup.rg.bam \\
			-O {config['output_dir']}/alignments/{config['FUO_sample_name']}.bwa.markdup.rg.bqsr.bam \\
			--bqsr-recal-file {config['output_dir']}/alignments/{config['FUO_sample_name']}.bwa.markdup.rg.bam.table"""

rule MU0_ApplyBQSR:
	input:	
		f"{config['output_dir']}/alignments/{config['MU0_sample_name']}.bwa.markdup.rg.bam.table"
	output:
		f"{config['output_dir']}/alignments/{config['MU0_sample_name']}.bwa.markdup.rg.bqsr.bam"
	shell:
		f"""{gatk} ApplyBQSRSpark -R {config['reference_panel_path']} \\
			-I {config['output_dir']}/alignments/{config['MU0_sample_name']}.bwa.markdup.rg.bam \\
			-O {config['output_dir']}/alignments/{config['MU0_sample_name']}.bwa.markdup.rg.bqsr.bam \\
			--bqsr-recal-file {config['output_dir']}/alignments/{config['MU0_sample_name']}.bwa.markdup.rg.bam.table"""

rule Child_ApplyBQSR:
	input:
		f"{config['output_dir']}/alignments/{config['CHILD_sample_name']}.bwa.markdup.rg.bam.table"
	output:
		f"{config['output_dir']}/alignments/{config['CHILD_sample_name']}.bwa.markdup.rg.bqsr.bam"
	shell:
		f"""{gatk} ApplyBQSRSpark -R {config['reference_panel_path']} \\
			-I {config['output_dir']}/alignments/{config['CHILD_sample_name']}.bwa.markdup.rg.bam \\
			-O {config['output_dir']}/alignments/{config['CHILD_sample_name']}.bwa.markdup.rg.bqsr.bam \\
			--bqsr-recal-file {config['output_dir']}/alignments/{config['CHILD_sample_name']}.bwa.markdup.rg.bam.table"""
# ====================== Index BQSR BAM ======================
rule FUO_ApplyBQSR_index_bam:
	input:
		f"{config['output_dir']}/alignments/{config['FUO_sample_name']}.bwa.markdup.rg.bqsr.bam"
	output:
		f"{config['output_dir']}/alignments/{config['FUO_sample_name']}.bwa.markdup.rg.bqsr.bam.bai"
	conda:
		"./first_step_mamba.yml"
	shell:
		f"samtools index {config['output_dir']}/alignments/{config['FUO_sample_name']}.bwa.markdup.rg.bqsr.bam" # gatk build bam index will create the index file for the bam file

rule MU0_ApplyBQSR_index_bam:
	input:
		f"{config['output_dir']}/alignments/{config['MU0_sample_name']}.bwa.markdup.rg.bqsr.bam"
	output:
		f"{config['output_dir']}/alignments/{config['MU0_sample_name']}.bwa.markdup.rg.bqsr.bam.bai"
	conda:
		"./first_step_mamba.yml"
	shell:
		f"samtools index {config['output_dir']}/alignments/{config['MU0_sample_name']}.bwa.markdup.rg.bqsr.bam" # gatk build bam index will create the index file for the bam file

rule Child_ApplyBQSR_index_bam:
	input:
		f"{config['output_dir']}/alignments/{config['CHILD_sample_name']}.bwa.markdup.rg.bqsr.bam"
	output:
		f"{config['output_dir']}/alignments/{config['CHILD_sample_name']}.bwa.markdup.rg.bqsr.bam.bai"
	conda:
		"./first_step_mamba.yml"
	shell:
		f"samtools index {config['output_dir']}/alignments/{config['CHILD_sample_name']}.bwa.markdup.rg.bqsr.bam" # gatk build bam index will create the index file for the bam file

# ====================== AnalyzeCovariates ======================
rule FUO_BaseRecalibrator2:
	input:
		f"{config['output_dir']}/alignments/{config['FUO_sample_name']}.bwa.markdup.rg.bqsr.bam.bai",
		f"{config['output_dir']}/alignments/{config['FUO_sample_name']}.bwa.markdup.rg.bqsr.bam",
		f"{config['reference_panel_path']}",
		f"{config['known_sites']}",
		f"{config['known_sites']}.tbi"
	output:
		f"{config['output_dir']}/alignments/{config['FUO_sample_name']}.bwa.markdup.rg.bqsr.bam.table"
	shell:
		f"""{gatk} BaseRecalibrator \\
			-R {config['reference_panel_path']} \\
			-I {config['output_dir']}/alignments/{config['FUO_sample_name']}.bwa.markdup.rg.bqsr.bam \\
			-O {config['output_dir']}/alignments/{config['FUO_sample_name']}.bwa.markdup.rg.bqsr.bam.table \\
			--known-sites {config['known_sites']}"""

rule FU0_AnalyzeCovariates:
	input:
		f"{config['output_dir']}/alignments/{config['FUO_sample_name']}.bwa.markdup.rg.bqsr.bam.table",
		f"{config['output_dir']}/alignments/{config['FUO_sample_name']}.bwa.markdup.rg.bam.table"
	output:
		f"{config['output_dir']}/alignments/{config['FUO_sample_name']}.bwa.markdup.rg.bqsr.bam.table.pdf"
	conda:
		"./gatk_R_plot.yml"
	shell:
		f"""{gatk} AnalyzeCovariates \\
			-before {config['output_dir']}/alignments/{config['FUO_sample_name']}.bwa.markdup.rg.bam.table \\
			-after {config['output_dir']}/alignments/{config['FUO_sample_name']}.bwa.markdup.rg.bqsr.bam.table \\
			-plots {config['output_dir']}/alignments/{config['FUO_sample_name']}.bwa.markdup.rg.bqsr.bam.table.pdf"""

rule MU0_BaseRecalibrator2:
	input:
		f"{config['output_dir']}/alignments/{config['MU0_sample_name']}.bwa.markdup.rg.bqsr.bam.bai",
		f"{config['output_dir']}/alignments/{config['MU0_sample_name']}.bwa.markdup.rg.bqsr.bam",
		f"{config['reference_panel_path']}",
		f"{config['known_sites']}",
		f"{config['known_sites']}.tbi"
	output:
		f"{config['output_dir']}/alignments/{config['MU0_sample_name']}.bwa.markdup.rg.bqsr.bam.table"
	shell:
		f"""{gatk} BaseRecalibrator \\
			-R {config['reference_panel_path']} \\
			-I {config['output_dir']}/alignments/{config['MU0_sample_name']}.bwa.markdup.rg.bqsr.bam \\
			-O {config['output_dir']}/alignments/{config['MU0_sample_name']}.bwa.markdup.rg.bqsr.bam.table \\
			--known-sites {config['known_sites']}"""

rule MU0_AnalyzeCovariates:
	input:
		f"{config['output_dir']}/alignments/{config['MU0_sample_name']}.bwa.markdup.rg.bqsr.bam.table",
		f"{config['output_dir']}/alignments/{config['MU0_sample_name']}.bwa.markdup.rg.bam.table"
	output:
		f"{config['output_dir']}/alignments/{config['MU0_sample_name']}.bwa.markdup.rg.bqsr.bam.table.pdf"
	conda:
		"./gatk_R_plot.yml"
	shell:
		f"""{gatk} AnalyzeCovariates \\
			-before {config['output_dir']}/alignments/{config['MU0_sample_name']}.bwa.markdup.rg.bam.table \\
			-after {config['output_dir']}/alignments/{config['MU0_sample_name']}.bwa.markdup.rg.bqsr.bam.table \\
			-plots {config['output_dir']}/alignments/{config['MU0_sample_name']}.bwa.markdup.rg.bqsr.bam.table.pdf"""

rule Child_BaseRecalibrator2:
	input:
		f"{config['output_dir']}/alignments/{config['CHILD_sample_name']}.bwa.markdup.rg.bqsr.bam.bai",
		f"{config['output_dir']}/alignments/{config['CHILD_sample_name']}.bwa.markdup.rg.bqsr.bam",
		f"{config['reference_panel_path']}",
		f"{config['known_sites']}",
		f"{config['known_sites']}.tbi"
	output:
		f"{config['output_dir']}/alignments/{config['CHILD_sample_name']}.bwa.markdup.rg.bqsr.bam.table"
	shell:
		f"""{gatk} BaseRecalibrator \\
			-R {config['reference_panel_path']} \\
			-I {config['output_dir']}/alignments/{config['CHILD_sample_name']}.bwa.markdup.rg.bqsr.bam \\
			-O {config['output_dir']}/alignments/{config['CHILD_sample_name']}.bwa.markdup.rg.bqsr.bam.table \\
			--known-sites {config['known_sites']}"""

rule Child_AnalyzeCovariates:
	input:
		f"{config['output_dir']}/alignments/{config['CHILD_sample_name']}.bwa.markdup.rg.bqsr.bam.table",
		f"{config['output_dir']}/alignments/{config['CHILD_sample_name']}.bwa.markdup.rg.bam.table"
	output:
		f"{config['output_dir']}/alignments/{config['CHILD_sample_name']}.bwa.markdup.rg.bqsr.bam.table.pdf"
	conda:
		"./gatk_R_plot.yml"
	shell:
		f"""{gatk} AnalyzeCovariates \\
			-before {config['output_dir']}/alignments/{config['CHILD_sample_name']}.bwa.markdup.rg.bam.table \\
			-after {config['output_dir']}/alignments/{config['CHILD_sample_name']}.bwa.markdup.rg.bqsr.bam.table \\
			-plots {config['output_dir']}/alignments/{config['CHILD_sample_name']}.bwa.markdup.rg.bqsr.bam.table.pdf"""

# ====================== HaplotypeCaller ======================
rule FU0_HaplotypeCaller:
	input:
		f"{config['output_dir']}/alignments/{config['FUO_sample_name']}.bwa.markdup.rg.bqsr.bam",
		f"{config['output_dir']}/alignments/{config['FUO_sample_name']}.bwa.markdup.rg.bqsr.bam.bai",
		f"{config['reference_panel_path']}"
	output:
		f"{config['output_dir']}/variants/{config['FUO_sample_name']}.g.vcf.gz"
	shell:
		f"""{gatk} HaplotypeCaller \\
			-R {config['reference_panel_path']} \\
			--emit-ref-confidence GVCF \\
			-I {config['output_dir']}/alignments/{config['FUO_sample_name']}.bwa.markdup.rg.bqsr.bam \\
			-O {config['output_dir']}/variants/{config['FUO_sample_name']}.g.vcf.gz \\
			-D {config['known_sites']} \\
			-ip 50 \\
			-L {config['bed_path']}
			"""

rule MU0_HaplotypeCaller:
	input:
		f"{config['output_dir']}/alignments/{config['MU0_sample_name']}.bwa.markdup.rg.bqsr.bam",
		f"{config['output_dir']}/alignments/{config['MU0_sample_name']}.bwa.markdup.rg.bqsr.bam.bai",
		f"{config['reference_panel_path']}"
	output:
		f"{config['output_dir']}/variants/{config['MU0_sample_name']}.g.vcf.gz"
	shell:
		f"""{gatk} HaplotypeCaller \\
			-R {config['reference_panel_path']} \\
			--emit-ref-confidence GVCF \\
			-I {config['output_dir']}/alignments/{config['MU0_sample_name']}.bwa.markdup.rg.bqsr.bam \\
			-O {config['output_dir']}/variants/{config['MU0_sample_name']}.g.vcf.gz \\
			-D {config['known_sites']} \\
			-ip 50 \\
			-L {config['bed_path']}
			"""

rule Child_HaplotypeCaller:
	input:
		f"{config['output_dir']}/alignments/{config['CHILD_sample_name']}.bwa.markdup.rg.bqsr.bam",
		f"{config['output_dir']}/alignments/{config['CHILD_sample_name']}.bwa.markdup.rg.bqsr.bam.bai",
		f"{config['reference_panel_path']}"
	output:
		f"{config['output_dir']}/variants/{config['CHILD_sample_name']}.g.vcf.gz"
	shell:
		f"""{gatk} HaplotypeCaller \\
			-R {config['reference_panel_path']} \\
			--emit-ref-confidence GVCF \\
			-I {config['output_dir']}/alignments/{config['CHILD_sample_name']}.bwa.markdup.rg.bqsr.bam \\
			-O {config['output_dir']}/variants/{config['CHILD_sample_name']}.g.vcf.gz \\
			-D {config['known_sites']} \\
			-ip 50 \\
			-L {config['bed_path']}
			"""

# ====================== GenotypeGVCFs ======================
rule FU0_GenotypeGVCFs:
	input:
		f"{config['output_dir']}/variants/{config['FUO_sample_name']}.g.vcf.gz",
		f"{config['reference_panel_path']}"
	output:
		f"{config['output_dir']}/variants/{config['FUO_sample_name']}.vcf.gz"
	shell:
		f"""{gatk} GenotypeGVCFs \\
			-R {config['reference_panel_path']} \\
			-V {config['output_dir']}/variants/{config['FUO_sample_name']}.g.vcf.gz \\
			-O {config['output_dir']}/variants/{config['FUO_sample_name']}.vcf.gz"""

rule MU0_GenotypeGVCFs:
	input:
		f"{config['output_dir']}/variants/{config['MU0_sample_name']}.g.vcf.gz",
		f"{config['reference_panel_path']}"
	output:
		f"{config['output_dir']}/variants/{config['MU0_sample_name']}.vcf.gz"
	shell:
		f"""{gatk} GenotypeGVCFs \\
			-R {config['reference_panel_path']} \\
			-V {config['output_dir']}/variants/{config['MU0_sample_name']}.g.vcf.gz \\
			-O {config['output_dir']}/variants/{config['MU0_sample_name']}.vcf.gz"""

rule Child_GenotypeGVCFs:
	input:
		f"{config['output_dir']}/variants/{config['CHILD_sample_name']}.g.vcf.gz",
		f"{config['reference_panel_path']}"
	output:
		f"{config['output_dir']}/variants/{config['CHILD_sample_name']}.vcf.gz"
	shell:
		f"""{gatk} GenotypeGVCFs \\
			-R {config['reference_panel_path']} \\
			-V {config['output_dir']}/variants/{config['CHILD_sample_name']}.g.vcf.gz \\
			-O {config['output_dir']}/variants/{config['CHILD_sample_name']}.vcf.gz"""

# ====================== CombineGVCFs ======================
rule CombineGVCFs:
	input:
		f"{config['output_dir']}/variants/{config['FUO_sample_name']}.g.vcf.gz",
		f"{config['output_dir']}/variants/{config['MU0_sample_name']}.g.vcf.gz",
		f"{config['output_dir']}/variants/{config['CHILD_sample_name']}.g.vcf.gz",
		f"{config['reference_panel_path']}",
		f"{config['known_sites']}"
	output:
		f"{config['output_dir']}/variants/Trio.g.vcf.gz"
	shell:
		f"""{gatk} CombineGVCFs \\
			-R {config['reference_panel_path']} \\
			-V {config['output_dir']}/variants/{config['FUO_sample_name']}.g.vcf.gz \\
			-V {config['output_dir']}/variants/{config['MU0_sample_name']}.g.vcf.gz \\
			-V {config['output_dir']}/variants/{config['CHILD_sample_name']}.g.vcf.gz \\
			-O {config['output_dir']}/variants/Trio.g.vcf.gz \\
			-ip 50 \\
			-L {config['bed_path']} \\
			-D {config['known_sites']}"""

# ====================== Trio GenotypeGVCFs ======================
rule Trio_GenotypeGVCFs:
	input:
		f"{config['output_dir']}/variants/Trio.g.vcf.gz",
		f"{config['reference_panel_path']}"
	output:
		f"{config['output_dir']}/variants/Trio.vcf.gz"
	shell:
		f"""{gatk} GenotypeGVCFs \\
			-R {config['reference_panel_path']} \\
			-V {config['output_dir']}/variants/Trio.g.vcf.gz \\
			-O {config['output_dir']}/variants/Trio.vcf.gz"""

# ====================== VariantFiltration ======================
rule SNP_SelectVariants:
	input:
		f"{config['output_dir']}/variants/Trio.vcf.gz",
		f"{config['reference_panel_path']}"
	output:
		f"{config['output_dir']}/variants/Trio.snp.vcf"
	shell:
		f"""{gatk} SelectVariants \\
			-R {config['reference_panel_path']} \\
			-V {config['output_dir']}/variants/Trio.vcf.gz \\
			--select-type-to-include SNP \\
			-O {config['output_dir']}/variants/Trio.snp.vcf
		"""

rule Indel_SelectVariants:
	input:
		f"{config['output_dir']}/variants/Trio.vcf.gz",
		f"{config['reference_panel_path']}"
	output:
		f"{config['output_dir']}/variants/Trio.indel.vcf"
	shell:
		f"""{gatk} SelectVariants \\
			-R {config['reference_panel_path']} \\
			-V {config['output_dir']}/variants/Trio.vcf.gz \\
			--select-type-to-include INDEL \\
			-O {config['output_dir']}/variants/Trio.indel.vcf
		"""

# ====================== VariantFiltration ======================
rule SNP_VariantFiltration:
	input:
		f"{config['output_dir']}/variants/Trio.snp.vcf"
	output:
		f"{config['output_dir']}/variants/Trio.snp.pass.vcf"
	shell:
		f"""{gatk} VariantFiltration
			-R {config['reference_panel_path']} \\
			-V {config['output_dir']}/variants/Trio.snp.vcf \\
			--filter-expression 'QUAL<30.0' --filter-name 'LOW_QUAL' \\
			--filter-expression 'QD<2.0' --filter-name 'LOW_QD' \\
			--filter-expression 'FS>60.0' --filter-name 'HIGH_FS' \\
			--filter-expression 'MQ<40.0' --filter-name 'LOW_MQ' \\
			--filter-expression 'MQRankSum<-12.5' --filter-name 'LOW_MQRS' \\
			--filter-expression 'ReadPosRankSum<-8.0' --filter-name 'LOW_RPRS' \\
			--filter-expression 'SOR>3.0' --filter-name 'HIGH_SOR' \\
			-O {config['output_dir']}/variants/Trio.snp.pass.vcf
		"""

rule Indel_VariantFiltration: # --filter-expression 'QUAL<30.0' --filter-name 'LOW_QUAL' --filter-expression  'QD<2.0' --filter-name 'LOW_QD' --filter-expression 'FS>200.0' --filter-name 'HIGH_FS' --filter-expression 'ReadPosRankSum<-20.0' --filter-name 'LOW_RPRS' --filter-expression 'SOR>10.0' --filter-name 'HIGH_SOR' 
	input:
		f"{config['output_dir']}/variants/Trio.indel.vcf"
	output:
		f"{config['output_dir']}/variants/Trio.indel.pass.vcf"
	shell:
		f"""{gatk} VariantFiltration
			-R {config['reference_panel_path']} \\
			-V {config['output_dir']}/variants/Trio.indel.vcf \\
			--filter-expression 'QUAL<30.0' --filter-name 'LOW_QUAL' \\
			--filter-expression 'QD<2.0' --filter-name 'LOW_QD' \\
			--filter-expression 'FS>200.0' --filter-name 'HIGH_FS' \\
			--filter-expression 'ReadPosRankSum<-20.0' --filter-name 'LOW_RPRS' \\
			--filter-expression 'SOR>10.0' --filter-name 'HIGH_SOR' \\
			-O {config['output_dir']}/variants/Trio.indel.pass.vcf
		"""

# ====================== Get filtered ======================
rule get_filtered_snp_VariantFiltration: #	"perl -ne 'chomp;if(\$_=~/^#/ || \$_ =~ /PASS/){print \"\$_\\n\"}' $Trio.snp.pass.vcf > $Trio.snp.filter.vcf\n",

	input:
		f"{config['output_dir']}/variants/Trio.snp.pass.vcf"
	output:
		f"{config['output_dir']}/variants/Trio.snp.filter.vcf"
	shell:
		f"""perl -ne 'chomp;if(\$_=~/^#/ || \$_ =~ /PASS/){print "\$_\\n"}' {config['output_dir']}/variants/Trio.snp.pass.vcf > {config['output_dir']}/variants/Trio.snp.filter.vcf"""

rule get_filtered_indel_VariantFiltration: #	"perl -ne 'chomp;if(\$_=~/^#/ || \$_ =~ /PASS/){print \"\$_\\n\"}' $Trio.indel.pass.vcf > $Trio.indel.filter.vcf\n"
	input:
		f"{config['output_dir']}/variants/Trio.indel.pass.vcf"
	output:
		f"{config['output_dir']}/variants/Trio.indel.filter.vcf"
	shell:
		f"""perl -ne 'chomp;if(\$_=~/^#/ || \$_ =~ /PASS/){print "\$_\\n"}' {config['output_dir']}/variants/Trio.indel.pass.vcf > {config['output_dir']}/variants/Trio.indel.filter.vcf"""

# ====================== FilterVar ======================
rule SNP_FilterVar:
	input:
		f"{config['output_dir']}/variants/Trio.snp.filter.vcf"
	output:
		f"{config['output_dir']}/variants/Trio.snp.annovar"
	shell:
		f"""perl {config['FilterVar_path']} \\
			-in {config['output_dir']}/variants/Trio.snp.filter.vcf \\
			-out {config['output_dir']}/variants/Trio.snp.annovar \\
			-minhet 0.20 --wildsample -qual 30 -dp 20 -adp 10 -gq 30 --rough
			"""

rule Indel_FilterVar:
	input:
		f"{config['output_dir']}/variants/Trio.indel.filter.vcf"
	output:
		f"{config['output_dir']}/variants/Trio.indel.annovar"
	shell:
		f"""perl {config['FilterVar_path']} \\
			-in {config['output_dir']}/variants/Trio.indel.filter.vcf \\
			-out {config['output_dir']}/variants/Trio.indel.annovar \\
			-minhet 0.30 --wildsample -qual 30 -dp 20 -adp 10 -gq 30 --rough
		""" # -minhet different from SNP_FilterVar

# ====================== Combine indel and snp ======================
rule Combine_snp_indel_annovar:
	input:
		f"{config['output_dir']}/variants/Trio.snp.annovar",
		f"{config['output_dir']}/variants/Trio.indel.annovar"
	output:
		f"{config['output_dir']}/variants/Trio.annovar"
	shell:
		f"""cat {config['output_dir']}/variants/Trio.snp.annovar {config['output_dir']}/variants/Trio.indel.annovar > {config['output_dir']}/variants/Trio.annovar"""

# ====================== Find denovo mutation in Trio ======================
rule DenovoCNN:# -w=$dir[2]/$Trio -cv=$dir[2]/$Proband/$Proband.vcf.gz -fv=$dir[2]/$Father/$Father.vcf.gz -mv=$dir[2]/$Mother/$Mother.vcf.gz -cb=$dir[1]/$Proband/$Proband.sort.marked.bam -fb=$dir[1]/$Father/$Father.sort.marked.bam -mb=$dir[1]/$Mother/$Mother.sort.marked.bam -g=$ref_fa -sm=$snp_model -im=$ins_model -dm=$del_model -o=$Proband",".DNMs_predictions.csv

	input:
		f"{config['output_dir']}/variants/{config['FUO_sample_name']}.vcf.gz",
		f"{config['output_dir']}/variants/{config['MU0_sample_name']}.vcf.gz",
		f"{config['output_dir']}/variants/{config['CHILD_sample_name']}.vcf.gz",
		f"{config['DenovoCNN_path']}"
	output:
		f"{config['output_dir']}/variants/{config['CHILD_sample_name']}.DNMs_predictions.csv"
	conda:
		f"{config['DenovoCNN_env']}"
	shell:
		f"""bash {config['DenovoCNN_path']} \\
			-w={config['output_dir']}/variants/{config['CHILD_sample_name']} \\
			-cv={config['output_dir']}/variants/{config['CHILD_sample_name']}.vcf.gz \\
			-fv={config['output_dir']}/variants/{config['FUO_sample_name']}.vcf.gz \\
			-mv={config['output_dir']}/variants/{config['MU0_sample_name']}.vcf.gz \\
			-cb={config['output_dir']}/alignments/{config['CHILD_sample_name']}.bwa.markdup.rg.bam \\
			-fb={config['output_dir']}/alignments/{config['FUO_sample_name']}.bwa.markdup.rg.bam \\
			-mb={config['output_dir']}/alignments/{config['MU0_sample_name']}.bwa.markdup.rg.bam \\
			-g={config['reference_panel_path']} \\
			-sm={config['snp_model']} \\
			-im={config['ins_model']} \\
			-dm={config['del_model']} \\
			-o={config['output_dir']}/variants/{config['CHILD_sample_name']}.DNMs_predictions.csv
			"""

rule Filter_DenovoCNN: # FilterDNMs -in $Proband",".DNMs_predictions.csv -out $Proband",".preDNMs.txt -p $opt{DNMs_p} -c $opt{DNMs_c}
	input:
		f"{config['output_dir']}/variants/{config['CHILD_sample_name']}.DNMs_predictions.csv"
	output:
		f"{config['output_dir']}/variants/{config['CHILD_sample_name']}.preDNMs.txt"
	shell:
		f"""perl {config['FilterDNMs_path']}
			-in {config['output_dir']}/variants/{config['CHILD_sample_name']}.DNMs_predictions.csv \\
			-out {config['output_dir']}/variants/{config['CHILD_sample_name']}.preDNMs.txt \\
			-p {config['DNMs_p']}
			-c {config['DNMs_c']}
		"""

# ====================== ExtremeVar ======================
rule ExtremeVar: # ExtremeVar -in $dir[2]/$Trio/$Trio.annovar -out $Trio.initial -Psoft $opt{Psoft} -maf $opt{MAFS} -reference $opt{ref} --extreme --extreme_all --remove -database $db_path  $db -TrioID $opt{TrioID}\n
	input:
		f"{config['ExtremeVar_PATH']}",
		f"{config['output_dir']}/variants/Trio.annovar"
	output:
		f"{config['output_dir']}/variants/Trio.initial"
	shell:
		f"""perl {config['ExtremeVar_PATH']} \\
		-in {config['output_dir']}/variants/Trio.annovar \\
		-out {config['output_dir']}/variants/Trio.initial \\
		-Psoft {config['Psoft']} \\
		-maf {config['MAFS']} \\
		-reference {config['reference_panel_name']} \\
		--extreme \\
		--extreme_all \\
		--remove \\
		-database {config['database_path']}  {db}
		-TrioID {Trio_ID}
		"""

# ====================== ExtremeVar2 ======================
rule ExtremeVar2:
	input:
		f"{config['ExtremeVar2_PATH']}",
		f"{config['output_dir']}/variants/Trio.initial.extreme.xls",
		f"{config['output_dir']}/variants/Trio.initial.hg38_multianno.txt"
	output:
		f"{config['output_dir']}/variants/Trio.initial.tmp.xls"
	shell:
		f"""perl {config['ExtremeVar2_PATH']} \\
		-out {config['output_dir']}/variants/Trio.initial \\
		--extreme \\
		--extreme_all \\
		-MPsoft {config['MPsoft2']} \\
		-reference {config['reference_panel_name']}
		"""

rule awk_select: # 		"awk -F \"\\t\" 'BEGIN{IGNORECASE=1} NR==1 {print \$0} NR>1 {if(\$1~/Y/ || \$105~/pathogenic/ || \$105~/drug_response/ || (\$107~/DM/ && \$105 !~/benign/)) print \$0}' $Trio.initial.tmp.xls | cut -f 1-10,12-26,30-33,37,41,45,51-53,55-60,64,69,79,85,88,93-94,97,100-111,115 > $Trio.initial.xls\n",
	input:
		f"{config['output_dir']}/variants/Trio.initial.tmp.xls"
	output:
		f"{config['output_dir']}/variants/Trio.initial.xls"
	shell:
		r"""awk -F \"\\t\" 'BEGIN{IGNORECASE=1} NR==1 {print \$0} NR>1 {if(\$1~/Y/ || \$105~/pathogenic/ || \$105~/drug_response/ || (\$107~/DM/ && \$105 !~/benign/)) print \$0}' {config['output_dir']}/variants/Trio.initial.tmp.xls | cut -f 1-10,12-26,30-33,37,41,45,51-53,55-60,64,69,79,85,88,93-94,97,100-111,115 > {config['output_dir']}/variants/Trio.initial.xls"""

# ====================== Process Trio ======================
rule Process_Trios: # $process_trio -in $Trio.initial.xls -TP $transcript -Inheritance $gm -HPO $hpo -TrioID $opt{TrioID} -out $Trio.extreme.xls
	input:
		f"{config['output_dir']}/variants/Trio.initial.xls"
	output:
		f"{config['output_dir']}/variants/Trio.extreme.xls"
	shell:
		f"""perl {config['process_trio_path']} \\
			-in {config['output_dir']}/variants/Trio.initial.xls \\
			-TP {config['TRANSCRIPT']} \\
			-Inheritance {config['GENEMODEL']} \\
			-HPO {config['HPO']} \\
			-TrioID {Trio_ID} \\
			-out {config['output_dir']}/variants/Trio.extreme.xls
		"""

# ====================== Sample ACMG ======================
# line 243
rule ACMG: # $acmg -in $Trio.extreme.xls -out $Trio.prefinal.xls 
	input:
		f"{config['ACMG_path']}",
		f"{config['output_dir']}/variants/Trio.extreme.xls"
	output:
		f"{config['output_dir']}/variants/Trio.prefinal.xls"
	shell:
		f"""perl {config['ACMG_path']} \\
			-in {config['output_dir']}/variants/Trio.extreme.xls \\
			-out {config['output_dir']}/variants/Trio.prefinal.xls
		"""
# ====================== Get trio ======================
rule get_Trio: # $getTrio $dir[2]/$Trio/$Proband",".preDNMs.txt $Trio.initial.tmp.xls $Trio.initial.Trios.extreme.xls | cut -f 1-11,13-27,31-34,38,42,46,52-54,56-61,65,70,80,86,89,94-95,98,101-112,116  > $Trio.initial.Trios.xls
	input:
		f"{config['output_dir']}/variants/{config['CHILD_sample_name']}.preDNMs.txt",
		f"{config['output_dir']}/variants/Trio.initial.tmp.xls"
	output:
		f"{config['output_dir']}/variants/Trio/Trio.initial.Trios.xls"
	shell:
		f"""perl {config['getTrio_path']} \\
			{config['output_dir']}/variants/{config['CHILD_sample_name']}.preDNMs.txt \\
			{config['output_dir']}/variants/Trio.initial.tmp.xls \\
			{config['output_dir']}/variants/Trio.initial.Trios.extreme.xls | cut -f 1-11,13-27,31-34,38,42,46,52-54,56-61,65,70,80,86,89,94-95,98,101-112,116  > {config['output_dir']}/variants/Trio/Trio.initial.Trios.xls
		"""

# ====================== Another process Trios =====================
rule process_Trios: # $process_Trios -in $Trio.initial.Trios.xls -TP $transcript -Inheritance $gm -HPO $hpo -TrioID $opt{TrioID} -out $Trio.Trios.xls

	input:
		f"{config['output_dir']}/variants/Trio/Trio.initial.Trios.xls"
	output:
		f"{config['output_dir']}/variants/Trio/Trio.Trios.xls"
	shell:
		f"""perl {config['process_Trios_path']} \\
			-in {config['output_dir']}/variants/Trio/Trio.initial.Trios.xls \\
			-TP {config['TRANSCRIPT']} \\
			-Inheritance {config['GENEMODEL']} \\
			-HPO {config['HPO']} \\
			-TrioID {Trio_ID} \\
			-out {config['output_dir']}/variants/Trio/Trio.Trios.xls
		"""

# ====================== Sample ACMG ======================
rule ACMG2: # $acmg -in $Trio.Trios.xls -out $Trio.Trios.filter.xls
	input:
		f"{config['ACMG_path']}",
		f"{config['output_dir']}/variants/Trio/Trio.Trios.xls"
	output:
		f"{config['output_dir']}/variants/Trio.Trios.filter.xls"
	shell:
		f"""perl {config['ACMG_path']} \\
			-in {config['output_dir']}/variants/Trio/Trio.Trios.xls \\
			-out {config['output_dir']}/variants/Trio.Trios.filter.xls
		"""

# ====================== OUTPUT ALL ======================
# """if($opt{type}=~/ASD/i){
# 		print SH4 "$getTrioVar -in $Trio.prefinal.xls -list $panel_path -acmglist $acmg78 -DNMs $Trio.Trios.filter.xls -o1 $Trio.clinvar_HGMD.xls -o2 $Trio.loose.xls -o3 $Trio.strict.xls -o4 $Trio.final.xls\n";
# 	}else{
# 		print SH4 "$getTrioVar -in $Trio.prefinal.xls -acmglist $acmg78 -DNMs $Trio.Trios.filter.xls -o1 $Trio.clinvar_HGMD.xls -o2 $Trio.loose.xls -o3 $Trio.strict.xls $Trio.final.xls\n";
# 	}
# """
rule getTrioVar:
	input:
		f"{config['getTrioVar_path']}",
		f"{config['output_dir']}/variants/Trio.prefinal.xls"
		f"{config['output_dir']}/variants/Trio.Trios.filter.xls"
	output:
		f"{config['output_dir']}/variants/Trio.clinvar_HGMD.xls",
		f"{config['output_dir']}/variants/Trio.loose.xls",
		f"{config['output_dir']}/variants/Trio.strict.xls",
		f"{config['output_dir']}/variants/Trio.final.xls"
	shell:
		f"""if [{getTrioVar_flag}]
			then
				perl {config['getTrioVar_path']} \\
					-in {config['output_dir']}/variants/Trio.prefinal.xls \\
					-list {config['panel_path']} \\
					-acmglist {config['acmg78_path']} \\
					-DNMs {config['output_dir']}/variants/Trio.Trios.filter.xls \\
					-o1 {config['output_dir']}/variants/Trio.clinvar_HGMD.xls \\
					-o2 {config['output_dir']}/variants/Trio.loose.xls \\
					-o3 {config['output_dir']}/variants/Trio.strict.xls \\
					-o4 {config['output_dir']}/variants/Trio.final.xls
			else
				perl {config['getTrioVar_path']} \\
					-in {config['output_dir']}/variants/Trio.prefinal.xls \\
					-acmglist {config['acmg78_path']} \\
					-DNMs {config['output_dir']}/variants/Trio.Trios.filter.xls \\
					-o1 {config['output_dir']}/variants/Trio.clinvar_HGMD.xls \\
					-o2 {config['output_dir']}/variants/Trio.loose.xls \\
					-o3 {config['output_dir']}/variants/Trio.strict.xls \\
					{config['output_dir']}/variants/Trio.final.xls
			fi
		"""

# ====================== Finally get result ======================
rule get_result: # $result $dir[3]/$Trio.Trios.filter.xls $dir[3]/$Trio.strict.xls $dir[3]/$Trio.clinvar_HGMD.xls $dir[3]/$Trio.loose.xls $dir[3]/$Trio.final.xls $opt{outdir}/$Trio.stat.xls $Trio.result.xls
	input:
		f"{config['get_result_path']}"
		f"{config['output_dir']}/variants/Trio.Trios.filter.xls",
		f"{config['output_dir']}/variants/Trio.strict.xls",
		f"{config['output_dir']}/variants/Trio.clinvar_HGMD.xls",
		f"{config['output_dir']}/variants/Trio.loose.xls",
		f"{config['output_dir']}/variants/Trio.final.xls"
	output:
		f"{config['output_dir']}/Result.xls"
	shell:
		f"""perl {config['get_result_path']} \\
			{config['output_dir']}/variants/Trio.Trios.filter.xls \\
			{config['output_dir']}/variants/Trio.strict.xls \\
			{config['output_dir']}/variants/Trio.clinvar_HGMD.xls \\
			{config['output_dir']}/variants/Trio.loose.xls \\
			{config['output_dir']}/variants/Trio.final.xls \\
			{config['output_dir']}/Result.xls
		"""


# ====================== Prepare environment ====================== 
rule index_reference:
	input:
		f"{config['reference_panel_path']}"
	output:
		f"{config['reference_panel_path']}.bwt",
		f"{config['reference_panel_path']}.pac",
		f"{config['reference_panel_path']}.ann",
		f"{config['reference_panel_path']}.sa",
		f"{config['reference_panel_path']}.amb"
	conda:
		"./first_step_mamba.yml"
	shell:
		f"bwa index {config['reference_panel_path']}" # bwa index will create the index files with the fasta file
		
rule dict_index:
	input:
		f"{config['reference_panel_path']}"
	output:
		f"{ref_dict_path}"
	shell:
		f"{gatk} CreateSequenceDictionary REFERENCE={config['reference_panel_path']} OUTPUT={ref_dict_path}"

rule IndexFeatureFile: # index the known sites(dbsnp, vcf file, gz)
	input:
		f"{config['known_sites']}"
	output:
		f"{config['known_sites']}.tbi"
	shell:
		f"""{gatk} IndexFeatureFile \\
			-I {config['known_sites']}"""

rule fasta_faidx:
	input:
		f"{config['reference_panel_path']}"
	output:
		f"{config['reference_panel_path']}.fai"
	conda:
		"./first_step_mamba.yml"
	shell:
		f"samtools faidx {config['reference_panel_path']}"
	# samtools faidx will create the index file for the fasta file

#==============================================================================================================