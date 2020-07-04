import pandas as pd
import os

localrules: all
#shell.prefix("singularity exec -B $SCRATCH /gpfs/fs0/scratch/n/nicholsa/zyfniu/nfcore-sarek-2.6.img ") ### to use singularity shell for all processes


SCRATCH = "/gpfs/fs0/scratch/n/nicholsa/zyfniu"
INPUT = pd.read_table( SCRATCH + "/Sample/HPV32.tsv",names = ['Patient','Sex','n_vs_t','Sample','Lane','Fastq1','Fastq2'])
INPUT['Lane'] = INPUT.Lane.apply(str)
INPUT['Sample_Lane'] = INPUT.Sample + "_" + INPUT.Lane
SAMPLE =  INPUT['Sample'].unique()
TMP = "temp"
PAT = INPUT.Patient.drop_duplicates()
TUMOR = INPUT[(INPUT.n_vs_t == 1)].Sample.drop_duplicates()
SAMPLE_LANE = INPUT.Sample_Lane

REF_fasta = "$SCRATCH/igenomes_ref/Homo_sapiens_assembly38.fasta"
REF_dbsnp = "$SCRATCH/igenomes_ref/dbsnp_146.hg38.vcf.gz"
REF_known_indels = "$SCRATCH/igenomes_ref/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
REF_gnomAD = "$SCRATCH/igenomes_ref/gnomAD.r2.1.1.GRCh38.PASS.AC.AF.only.vcf.gz"
REF_pon = "$SCRATCH/igenomes_ref/1000g_pon.hg38.vcf.gz"
REF_exac_common = "/gpfs/fs0/scratch/n/nicholsa/zyfniu/igenomes_ref/somatic-hg38_small_exac_common_3.hg38.vcf.gz"

CHROMOSOMES = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
               'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19',
               'chr20', 'chr21', 'chr22', 'chrX', 'chrY']


OUTDIR = "AN_WGS/"


rule all:
	input:
		expand("QC/{sample_lane}/{sample_lane}_R1_fastqc.html", sample_lane = SAMPLE_LANE),
		expand("orphan/{sample}/Recal/{sample}.recal.bam", sample = SAMPLE),
		expand("results/mutect2/{patient}/filtered_{tumor}_vs_{patient}-N.vcf.idx",zip, patient = [x[:-2] for x in TUMOR],
		tumor = TUMOR) ### using zip https://endrebak.gitbooks.io/the-snakemake-book/chapters/expand/expand.html

def bwa_mem_fastq1(wildcards):
	return expand(INPUT[INPUT.Sample_Lane == wildcards.sample_lane].Fastq1)

def bwa_mem_fastq2(wildcards):
	return expand(INPUT[INPUT.Sample_Lane == wildcards.sample_lane].Fastq2)

rule fastqc:
	input:
        	fastq1 = bwa_mem_fastq1,
        	fastq2 = bwa_mem_fastq2
	output:
		o1 =  "QC/{sample_lane}/{sample_lane}_R1_fastqc.html",
		o2 =  "QC/{sample_lane}/{sample_lane}_R2_fastqc.html",
		o3 =  "QC/{sample_lane}/{sample_lane}_R1_fastqc.zip",
		o4 =  "QC/{sample_lane}/{sample_lane}_R2_fastqc.zip"
	threads: 10
	group: "align"
	shell:
        	"""
 		singularity exec -B $SCRATCH /gpfs/fs0/scratch/n/nicholsa/zyfniu/nfcore-sarek-2.6.img \
		fastqc -t {threads} {input} --outdir=QC/{wildcards.sample_lane}/
        mv QC/{wildcards.sample_lane}/*R1*fastqc.html QC/{wildcards.sample_lane}/{wildcards.sample_lane}_R1_fastqc.html
		mv QC/{wildcards.sample_lane}/*R2*fastqc.html QC/{wildcards.sample_lane}/{wildcards.sample_lane}_R2_fastqc.html
		mv QC/{wildcards.sample_lane}/*R1*fastqc.zip QC/{wildcards.sample_lane}/{wildcards.sample_lane}_R1_fastqc.zip
		mv QC/{wildcards.sample_lane}/*R2*fastqc.zip QC/{wildcards.sample_lane}/{wildcards.sample_lane}_R2_fastqc.zip
        	"""

def createRG(wildcards):
		return expand("@RG\\tID:{idRun}\\tPU:{idRun}\\tSM:{idSample}\\tLB:{idSample}\\tPL:illumina",
		idRun = INPUT[INPUT.Sample_Lane == wildcards.sample_lane].Lane,
		idSample = INPUT[INPUT.Sample_Lane == wildcards.sample_lane].Sample)

rule bwa_mem:
	input:
		fastq1 = bwa_mem_fastq1,
		fastq2 = bwa_mem_fastq2
		#fastq1= INPUT[INPUT.Sample_Lane == {sample_lane}].Fastq1,
		#fastq2= INPUT[INPUT.Sample_Lane == {sample_lane}].Fastq2
	output:
		temp("orphan/bwa/{sample_lane}.bwa.bam")
	threads: 70
	group: "align"
	params:
		readGroup = createRG
	shell:
#		INPUT.loc[INPUT['Sample_Lane']=={wildcards.sample_lane},'Bam'] = {output}
		"""
		singularity exec -B $SCRATCH /gpfs/fs0/scratch/n/nicholsa/zyfniu/nfcore-sarek-2.6.img \
		bwa mem -K 100000000 -R \"{params.readGroup}\" -B 3 -t {threads} -M {REF_fasta} \
    		{input.fastq1} {input.fastq2} | singularity exec -B $SCRATCH /gpfs/fs0/scratch/n/nicholsa/zyfniu/nfcore-sarek-2.6.img \
		samtools sort --threads {threads} -m 2G - > {output}
		"""


def get_bams_to_merge(wildcards):
	SAMPLE_LANE = INPUT[INPUT.Sample == wildcards.sample].Sample_Lane
	return expand("orphan/bwa/{sample_lane}.bwa.bam", sample_lane = SAMPLE_LANE)

rule merge_bam_mapped_and_index:
	input:
		get_bams_to_merge
	output:
		temp("orphan/bwa/{sample}.merged.bam")
	threads: 10
	group: "merge_markduplicate"
	run:
		if len(input) == 1:
			shell("mv {input} {output}")
			shell("singularity exec -B $SCRATCH /gpfs/fs0/scratch/n/nicholsa/zyfniu/nfcore-sarek-2.6.img samtools index {output}")
		elif len(input) > 1:
			#shell("mv {input} {output}")
			shell("singularity exec -B $SCRATCH /gpfs/fs0/scratch/n/nicholsa/zyfniu/nfcore-sarek-2.6.img samtools merge --threads {threads} - {input} | tee {output} | singularity exec -B $SCRATCH /gpfs/fs0/scratch/n/nicholsa/zyfniu/nfcore-sarek-2.6.img samtools index -@ {threads} - ")

rule markdup:
	input:
		"orphan/bwa/{sample}.merged.bam"
	output:
		bam = "orphan/MD/{sample}.md.bam",
		metric = "orphan/MD/{sample}.bam.metric"
	threads: 20
	group: "merge_markduplicate"
	shell:
		"""
		singularity exec -B $SCRATCH /gpfs/fs0/scratch/n/nicholsa/zyfniu/nfcore-sarek-2.6.img \
		gatk --java-options "-Xmx16G" \
		MarkDuplicates \
		--MAX_RECORDS_IN_RAM 50000 \
		--INPUT {input} \
		--METRICS_FILE {output.metric} \
		--TMP_DIR {TMP} \
		--ASSUME_SORT_ORDER coordinate \
		--CREATE_INDEX true \
		--OUTPUT {output.bam}
		"""

rule recalibrator:
	input:
		"orphan/MD/{sample}.md.bam"
	output:
		temp("orphan/{sample}/Recal/{sample}.{chr}.recal.table")
	threads: 2
	group: "recalibrator"
	shell:
		"""
		singularity exec -B $SCRATCH /gpfs/fs0/scratch/n/nicholsa/zyfniu/nfcore-sarek-2.6.img \
		gatk --java-options "-Xmx8g" \
		BaseRecalibrator \
		-I {input} \
		-O {output} \
		--tmp-dir /tmp \
		-R {REF_fasta} \
		-L {wildcards.chr} \
		--known-sites {REF_dbsnp} \
		--known-sites {REF_known_indels} \
		--verbosity INFO
		"""

def recal_tbl_to_gather(wildcards):
	return expand("orphan/" + wildcards.sample + "/Recal/" + wildcards.sample + ".{chr}.recal.table",chr = CHROMOSOMES)

rule gather_recal_tbl:
	input:
		recal_tbl_to_gather
	output:
		"orphan/{sample}/Recal/{sample}.recal.table"
	threads: 2
	group: "recalibrator"
	run:
		command = 'singularity exec -B $SCRATCH /gpfs/fs0/scratch/n/nicholsa/zyfniu/nfcore-sarek-2.6.img gatk --java-options \"-Xmx8g\" GatherBQSRReports'
		for i in input:
			command = command + " -I " + i
		command = command + ' -O {output}'
		shell(command)

rule ApplyBQSR:
	input:
		bam = "orphan/MD/{sample}.md.bam",
		recal_table = "orphan/{sample}/Recal/{sample}.recal.table"
	output:
		temp("orphan/{sample}/Recal/{sample}.recal.{chr}.bam")
	group: "recalibrator"
	shell:
		"""
		singularity exec -B $SCRATCH /gpfs/fs0/scratch/n/nicholsa/zyfniu/nfcore-sarek-2.6.img \
		gatk --java-options "-Xmx8G" \
		ApplyBQSR \
		-R {REF_fasta} \
		--input {input.bam} \
		--output {output} \
		-L {wildcards.chr} \
		--bqsr-recal-file {input.recal_table}
		"""

def recal_bam_to_gather(wildcards):
	return expand("orphan/" + wildcards.sample + "/Recal/" + wildcards.sample + ".recal.{chr}.bam",chr = CHROMOSOMES)

rule Merge_Recal_Bam_and_index:
	input:
		recal_bam_to_gather
	output:
		bam = "orphan/{sample}/Recal/{sample}.recal.bam",
		index = "orphan/{sample}/Recal/{sample}.recal.bai"
	group: "recalibrator"
	threads: 20
	shell:
		"""
		singularity exec -B $SCRATCH /gpfs/fs0/scratch/n/nicholsa/zyfniu/nfcore-sarek-2.6.img \
		samtools merge --threads {threads} {output.bam} {input}
		singularity exec -B $SCRATCH /gpfs/fs0/scratch/n/nicholsa/zyfniu/nfcore-sarek-2.6.img \
		samtools index {output.bam} {output.index}
		"""

rule mutect2_multi:
	input:
		normal = "orphan/{patient}-N/Recal/{patient}-N.recal.bam",
		tumor = "orphan/{tumor}/Recal/{tumor}.recal.bam"
	#	recurrence = "{sample}R/Recal/{sample}R.recal.bam"
	output:
		vcf = temp("results/mutect2/{patient}/{tumor}_vs_{patient}-N.{chr}.vcf"),
		f12 = temp("results/mutect2/{patient}/{tumor}_vs_{patient}-N_f12.{chr}.tar.gz")
	threads: 2
	group: "mutect2"
	shell:
		"""
		singularity exec -B $SCRATCH /gpfs/fs0/scratch/n/nicholsa/zyfniu/nfcore-sarek-2.6.img \
		gatk --java-options "-Xmx8g" Mutect2 -R {REF} \
		-I {input.normal}  \
		-I {input.tumor} \
		-normal {wildcards.patient}N \
		-L {wildcards.chr} \
		--germline-resource {GERMLINE} \
		--panel-of-normals {PON} --f1r2-tar-gz {output.f12} \
		-O {output.vcf}
		"""

def concat_vcf(wildcards):
	return expand("results/mutect2/" + wildcards.patient + "/" + wildcards.tumor + "_vs_" + wildcards.patient + "-N.{chr}.vcf", chr = CHROMOSOMES)

rule merge_mutect2_vcf:
	input:
		concat_vcf
	output:
		"results/mutect2/{sample}/{tumor}_vs_{patient}-N.vcf"
	threads: 2
	group: "mutect2"
	shell:
		"""
		singularity exec -B $SCRATCH /gpfs/fs0/scratch/n/nicholsa/zyfniu/nfcore-sarek-2.6.img \
		vcf-concat {input} > {output}
		"""

def concat_vcf_stats(wildcards):
	return expand("results/mutect2/" + wildcards.patient + "/" + wildcards.tumor + "_vs_" + wildcards.patient + "-N.{chr}.vcf.stats", chr = CHROMOSOMES)

rule merge_stats:
	input:
		concat_vcf_stats
	output:
		"results/mutect2/{patient}/{tumor}_vs_{patient}-N_merged.stats"
	threads: 2
	group: "mutect2"
	run:
		import os, glob
		cmd = "singularity exec -B $SCRATCH /gpfs/fs0/scratch/n/nicholsa/zyfniu/nfcore-sarek-2.6.img gatk MergeMutectStats "
		for input_file in input:
			cmd = cmd + " --stats " + input_file
		cmd = cmd + " -O {output}"
		shell(cmd)

def concat_vcf_f12(wildcards):
	return expand("results/mutect2/" + wildcards.patient + "/" + wildcards.tumor + "_vs_" + wildcards.patient + "-N_f12.{chr}.tar.gz", chr = CHROMOSOMES)

rule gatk_LearnOrientationModel:
	input:
		concat_vcf_f12
	output:
		model = "results/mutect2/{patient}/{tumor}_read_orientation_model.tar.gz"
	group: "mutect2"
	run:
		import os, glob
		cmd = "singularity exec -B $SCRATCH /gpfs/fs0/scratch/n/nicholsa/zyfniu/nfcore-sarek-2.6.img gatk LearnReadOrientationModel "
		for input_file in input:
			cmd = cmd + " -I " + input_file
		cmd = cmd + " -O {output}"
		shell(cmd)

	#	gatk --java-options "-Xmx8g" LearnReadOrientationModel \
	#	-I {input} -O {output.model}
	#	"""

rule gatk_get_pileupsummaries_normal:
	input:
		##Downloaded from https://console.cloud.google.com/storage/browser/_details/gatk-best-practices/somatic-hg38/
		common_biallelic_vcf =  REF_exac_common,
		normal = "orphan/{patient}-N/Recal/{patient}-N.recal.bam"
	output:
		summary = "results/mutect2/{patient}/{patient}-N_pileupsummaries.table"
	threads: 2
	group: "mutect2"
	shell:
		"""
		singularity exec -B $SCRATCH /gpfs/fs0/scratch/n/nicholsa/zyfniu/nfcore-sarek-2.6.img \
		gatk --java-options "-Xmx8g" GetPileupSummaries \
		-I {input.normal} \
		-V {input.common_biallelic_vcf} \
		-L {input.common_biallelic_vcf} \
		-O {output.summary}
		"""

rule gatk_get_pileupsummaries_tumor:
	input:
		##Downloaded from https://console.cloud.google.com/storage/browser/_details/gatk-best-practices/somatic-hg38/
		common_biallelic_vcf =  REF_exac_common,
		primary = "orphan/{tumor}/Recal/{tumor}.recal.bam"
	output:
		summary = "results/mutect2/{patient}/{tumor}_pileupsummaries.table"
	threads: 2
	group: "mutect2"
	shell:
		"""
		singularity exec -B $SCRATCH /gpfs/fs0/scratch/n/nicholsa/zyfniu/nfcore-sarek-2.6.img \
		gatk --java-options "-Xmx8g" GetPileupSummaries \
		-I {input.primary} \
		-V {input.common_biallelic_vcf} \
		-L {input.common_biallelic_vcf} \
		-O {output.summary}
		"""

rule gatk_calcContam_primary:
	input:
		normal = "results/mutect2/{patient}/{patient}-N_pileupsummaries.table",
		primary = "results/mutect2/{patient}/{tumor}_pileupsummaries.table"
	output:
		contamTable = "results/mutect2/{patient}/{tumor}_calContam.table",
		segment = "results/mutect2/{patient}/{tumor}_tumor.segment"
	threads: 2
	group: "mutect2"
	shell:
		"""
		singularity exec -B $SCRATCH /gpfs/fs0/scratch/n/nicholsa/zyfniu/nfcore-sarek-2.6.img \
		gatk --java-options "-Xmx8g" CalculateContamination \
		-I {input.normal} \
		-matched {input.primary} \
		--tumor-segmentation {output.segment} \
		-O {output.contamTable}
		"""

rule index_unfiltered_vcf:
	input:
		vcf = "results/mutect2/{patient}/{tumor}_vs_{patient}-N.vcf",
	output:
		vcf_tbi = "results/mutect2/{patient}/{tumor}_vs_{patient}-N.vcf.idx"
	threads: 2
	group: "mutect2"
	shell:
		"""
		singularity exec -B $SCRATCH /gpfs/fs0/scratch/n/nicholsa/zyfniu/nfcore-sarek-2.6.img \
		gatk --java-options "-Xmx8g" IndexFeatureFile \
		-I {input.vcf} \
		--output {output.vcf_tbi}
		"""

rule gatk_filterMutect:
	input:
		#p_contamTable = "results/mutect2/{patient}/{tumor}_calContam.table",
		#p_segment = "results/mutect2/{patient}/{tumor}_tumor.segment",
		vcf_tbi = "results/mutect2/{patient}/{tumor}_vs_{patient}-N.vcf.idx",
		# r_contamTable = "results/mutect2/{patient}/recurrence_calContam.table",
		# r_segment = "results/mutect2/{patient}/recurrence_tumor.segment",
		vcf = "results/mutect2/{patient}/{tumor}_vs_{patient}-N.vcf",
		model = "results/mutect2/{patient}/{tumor}_read_orientation_model.tar.gz",
		stats = "results/mutect2/{patient}/{tumor}_vs_{patient}-N_merged.stats"
	output:
		filter_vcf = temp("results/mutect2/{patient}/filtered_{tumor}_vs_{patient}-N.{chr}.vcf")
	threads: 2
	group: "mutect2"
	shell:
		"""
		singularity exec -B $SCRATCH /gpfs/fs0/scratch/n/nicholsa/zyfniu/nfcore-sarek-2.6.img \
		gatk --java-options "-Xmx8g" FilterMutectCalls \
		--ob-priors {input.model} \
		-stats {input.stats} \
		-R {REF} \
		-V {input.vcf} \
		-L {wildcards.chr} \
		--output {output.filter_vcf}
		"""

def concat_vcf_filtered(wildcards):
	return	expand("results/mutect2/" + wildcards.patient + "/filtered_" + wildcards.tumor + "_vs_" + wildcards.patient + "--N.{chr}.vcf", chr = CHROMOSOMES)

rule merge_mutect2_vcf_filtered:
	input:
		concat_vcf_filtered
	output:
		"results/mutect2/{patient}/filtered_{tumor}_vs_{patient}-N.vcf"
	threads: 2
	group: "mutect2"
	shell:
		"singularity exec -B $SCRATCH /gpfs/fs0/scratch/n/nicholsa/zyfniu/nfcore-sarek-2.6.img vcf-concat {input} > {output}"


rule index_filtered_vcf:
	input:
		vcf = "results/mutect2/{patient}/filtered_{tumor}_vs_{patient}-N.vcf",
	output:
		vcf_tbi = "results/mutect2/{patient}/filtered_{tumor}_vs_{patient}-N.vcf.idx"
	threads: 2
	group: "mutect2"
	shell:
		"""
		singularity exec -B $SCRATCH /gpfs/fs0/scratch/n/nicholsa/zyfniu/nfcore-sarek-2.6.img \
		gatk --java-options "-Xmx8g" IndexFeatureFile \
		-I {input.vcf} \
		--output {output.vcf_tbi}
		"""

# rule annotate:
# 	input:
# 		vcf = "results/mutect2/{patient}/forced_{patient}R_vs_{patient}-N.{chr}.vcf"
# 	output:
#
# 	threads: 2
# 	group: "mutect2"
# 	shell:
# 		"""
# 		singularity exec -B $SCRATCH /gpfs/fs0/scratch/n/nicholsa/zyfniu/nfcore-sarek-2.6.img \
# 		snpEff -Xmx8g \
#         ${snpeffDb} \
#         -csvStats ${reducedVCF}_snpEff.csv \
#         -nodownload \
#         ${cache} \
#         -canon \
#         -v \
#         ${vcf} \
#         > ${reducedVCF}_snpEff.ann.vcf
# 		mv snpEff_summary.html ${reducedVCF}_snpEff.html
# 		bgzip < ${vcf} > ${vcf}.gz
# 		tabix ${vcf}.gz
# 		"""
