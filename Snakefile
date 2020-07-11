
import pandas as pd
import os

#localrules: all

###
### LOAD SAMPLES
###
SCRATCH = "/gpfs/fs0/scratch/n/nicholsa/zyfniu"
INPUT = pd.read_table("/scratch/n/nicholsa/zyfniu/AN_WGS/AN_WGS/Sample/Sample.tsv",names = ['Patient','Sex','n_vs_t','Sample','Lane','Fastq1','Fastq2'])
INPUT['Lane'] = INPUT.Lane.apply(str)
INPUT['Sample_Lane'] = INPUT.Sample + "_" + INPUT.Lane
SAMPLE =  INPUT['Sample'].unique()
TMP = "temp"
PAT = INPUT.Patient.drop_duplicates()
TUMOR = INPUT[(INPUT.n_vs_t == 1)].Sample.drop_duplicates()
SAMPLE_LANE = INPUT.Sample_Lane
print(INPUT)

###
###	REFERENCE FILE
###
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

###
### Final Results
###
rule all:
	input:
		###FASTQC + BWA Alignment
		expand("QC/{sample_lane}/{sample_lane}_R1_fastqc.html", sample_lane = SAMPLE_LANE),
		expand("/scratch/n/nicholsa/zyfniu/AN_WGS/orphan/{sample}/Recal/{sample}.recal.bam", sample = SAMPLE),
		### using zip https://endrebak.gitbooks.io/the-snakemake-book/chapters/expand/expand.html
		###BAMQC + samtools_stats
		expand("QC/{sample}/{sample}.samtools.stats.out",sample = SAMPLE),
		expand("QC/{sample}/bamQC/qualimapReport.html",sample = SAMPLE),
		####
		####VARIANT CALLING OUTPUTS
		####
		expand("/scratch/n/nicholsa/zyfniu/AN_WGS/results/mutect2/{tumor}_vs_{patient}-N/{tumor}_vs_{patient}-N_snpEff.ann.vcf.gz",zip, patient = [x[:-2] for x in TUMOR],tumor = TUMOR),
		#### Manta
		expand("/scratch/n/nicholsa/zyfniu/AN_WGS/results/Manta/{tumor}_vs_{patient}-N/Manta_snpeff_{tumor}_vs_{patient}-N.candidateSV.ann.vcf.gz",zip, patient = [x[:-2] for x in TUMOR],tumor = TUMOR),
		expand("/scratch/n/nicholsa/zyfniu/AN_WGS/results/Manta/{tumor}_vs_{patient}-N/Manta_snpeff_{tumor}_vs_{patient}-N.candidateSmallIndels.ann.vcf.gz",zip, patient = [x[:-2] for x in TUMOR],tumor = TUMOR),
		expand("/scratch/n/nicholsa/zyfniu/AN_WGS/results/Manta/{tumor}_vs_{patient}-N/Manta_snpeff_{tumor}_vs_{patient}-N.diploidSV.ann.vcf.gz",zip, patient = [x[:-2] for x in TUMOR],tumor = TUMOR),
		expand("/scratch/n/nicholsa/zyfniu/AN_WGS/results/Manta/{tumor}_vs_{patient}-N/Manta_snpeff_{tumor}_vs_{patient}-N.somaticSV.ann.vcf.gz",zip, patient = [x[:-2] for x in TUMOR],tumor = TUMOR),
		#### ASCAT
		expand("/scratch/n/nicholsa/zyfniu/AN_WGS/results/ASCAT/{tumor}_vs_{patient}-N/{tumor}_vs_{patient}-N.tumor.cnvs.txt",zip, tumor = TUMOR, patient = [x[:-2] for x in TUMOR])
		#### QC
		#expand("reports/{patient}_multiqc.html",patient = PAT)
	threads: 80
	shell:
		"singularity exec $SCRATCH/singularity_images/nfcore-sarek-2.6.img multiqc . -n reports/multiqc.html"
###
###	Step 1: BWA-MEM Alignment
###

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
	threads: 20
	group: "align"
	shell:
			"""
		singularity exec -B $SCRATCH/igenomes_ref,$SCRATCH/HPV_WGS/raw /gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/nfcore-sarek-2.6.img \
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
		"""
		singularity exec -B $SCRATCH/igenomes_ref,$SCRATCH/HPV_WGS/raw /gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/nfcore-sarek-2.6.img \
		bwa mem -K 100000000 -R \"{params.readGroup}\" -B 3 -t {threads} -M {REF_fasta} \
			{input.fastq1} {input.fastq2} | singularity exec -B $SCRATCH/igenomes_ref,$SCRATCH/HPV_WGS/raw /gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/nfcore-sarek-2.6.img \
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
			shell("singularity exec -B $SCRATCH/igenomes_ref,$SCRATCH/HPV_WGS/raw /gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/nfcore-sarek-2.6.img samtools index {output}")
		elif len(input) > 1:
			#shell("mv {input} {output}")
			shell("singularity exec -B $SCRATCH/igenomes_ref,$SCRATCH/HPV_WGS/raw /gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/nfcore-sarek-2.6.img samtools merge --threads {threads} - {input} | tee {output} | singularity exec -B $SCRATCH/igenomes_ref,$SCRATCH/HPV_WGS/raw /gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/nfcore-sarek-2.6.img samtools index -@ {threads} - ")

rule markdup:
	input:
		"orphan/bwa/{sample}.merged.bam"
	output:
		bam = temp("orphan/MD/{sample}.md.bam"),
		metric = "orphan/MD/{sample}.bam.metric"
	threads: 20
	group: "merge_markduplicate"
	shell:
		"""
		singularity exec -B $SCRATCH/igenomes_ref,$SCRATCH/HPV_WGS/raw /gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/gatk-4.1.8.img \
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
		singularity exec -B $SCRATCH/igenomes_ref,$SCRATCH/HPV_WGS/raw /gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/gatk-4.1.8.img \
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
		temp("orphan/{sample}/Recal/{sample}.recal.table")
	threads: 2
	group: "recalibrator"
	run:
		command = 'singularity exec -B $SCRATCH/igenomes_ref,$SCRATCH/HPV_WGS/raw /gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/nfcore-sarek-2.6.img gatk --java-options \"-Xmx8g\" GatherBQSRReports'
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
		singularity exec -B $SCRATCH/igenomes_ref,$SCRATCH/HPV_WGS/raw /gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/gatk-4.1.8.img \
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
		bam = "/scratch/n/nicholsa/zyfniu/AN_WGS/orphan/{sample}/Recal/{sample}.recal.bam",
		index = "/scratch/n/nicholsa/zyfniu/AN_WGS/orphan/{sample}/Recal/{sample}.recal.bai"
	group: "recalibrator"
	threads: 20
	shell:
		"""
		singularity exec -B $SCRATCH/igenomes_ref,$SCRATCH/HPV_WGS/raw /gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/nfcore-sarek-2.6.img \
		samtools merge --threads {threads} {output.bam} {input}
		singularity exec -B $SCRATCH/igenomes_ref,$SCRATCH/HPV_WGS/raw /gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/nfcore-sarek-2.6.img \
		samtools index {output.bam} {output.index}
		"""

###
### QC On BAM FILES
###
rule samtools_stats:
	input:
		bam = "/scratch/n/nicholsa/zyfniu/AN_WGS/orphan/{sample}/Recal/{sample}.recal.bam",
		index = "orphan/{sample}/Recal/{sample}.recal.bai"
	output:
			stats = "/scratch/n/nicholsa/zyfniu/AN_WGS/QC/{sample}/{sample}.samtools.stats.out"
	group: "qc"
	threads: 2
	shell:
		"""
		singularity exec /gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/nfcore-sarek-2.6.img \
		samtools stats {input.bam} > {output}
		"""

rule bamqc:
	input:
		bam = "/scratch/n/nicholsa/zyfniu/AN_WGS/orphan/{sample}/Recal/{sample}.recal.bam",
		index = "/scratch/n/nicholsa/zyfniu/AN_WGS/orphan/{sample}/Recal/{sample}.recal.bai"
	output:
			stats = "QC/{sample}/bamQC/qualimapReport.html"
	group: "qc"
	threads: 40
	shell:
		"""
		singularity exec /gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/nfcore-sarek-2.6.img \
		qualimap --java-mem-size=100G \
		bamqc \
		-bam {input.bam} \
		--paint-chromosome-limits \
		--genome-gc-distr HUMAN \
		-nt {threads} \
		-skip-duplicated \
		--skip-dup-mode 0 \
		-outdir QC/{wildcards.sample}/bamQC/ \
		-outformat HTML
		"""

###
### VARIANT CALLING
###

rule mutect2:
	input:
		normal = "/scratch/n/nicholsa/zyfniu/AN_WGS/orphan/{patient}-N/Recal/{patient}-N.recal.bam",
		tumor = "/scratch/n/nicholsa/zyfniu/AN_WGS/orphan/{tumor}/Recal/{tumor}.recal.bam"
	#	recurrence = "{sample}R/Recal/{sample}R.recal.bam"
	output:
		vcf = temp("/scratch/n/nicholsa/zyfniu/AN_WGS/results/mutect2/{tumor}_vs_{patient}-N/unfiltered_{tumor}_vs_{patient}-N.{chr}.vcf"),
		stats = temp("/scratch/n/nicholsa/zyfniu/AN_WGS/results/mutect2/{tumor}_vs_{patient}-N/unfiltered_{tumor}_vs_{patient}-N.{chr}.vcf.stats"),
		f12 = temp("/scratch/n/nicholsa/zyfniu/AN_WGS/results/mutect2/{tumor}_vs_{patient}-N/unfiltered_{tumor}_vs_{patient}-N_f12.{chr}.tar.gz")
	threads: 2
	group: "variantCalling"
	shell:
		"""
		singularity exec -B $SCRATCH/igenomes_ref,$SCRATCH/HPV_WGS/raw /gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/gatk-4.1.8.img \
		gatk --java-options "-Xmx8g" Mutect2 -R {REF_fasta} \
		-I {input.normal}  \
		-I {input.tumor} \
		-normal {wildcards.patient}-N \
		-L {wildcards.chr} \
		--germline-resource {REF_gnomAD} \
		--panel-of-normals {REF_pon} --f1r2-tar-gz {output.f12} \
		-O {output.vcf}
		"""

def concat_vcf(wildcards):
	return expand("/scratch/n/nicholsa/zyfniu/AN_WGS/results/mutect2/" + wildcards.tumor + "_vs_" + wildcards.patient + "-N" + "/unfiltered_" + wildcards.tumor + "_vs_" + wildcards.patient + "-N.{chr}.vcf", chr = CHROMOSOMES)

rule merge_mutect2_vcf:
	input:
		concat_vcf
	output:
		temp("/scratch/n/nicholsa/zyfniu/AN_WGS/results/mutect2/{tumor}_vs_{patient}-N/unfiltered_{tumor}_vs_{patient}-N.vcf")
	threads: 2
	group: "variantCalling"
	shell:
		"""
		singularity exec -B $SCRATCH/igenomes_ref,$SCRATCH/HPV_WGS/raw /gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/vcftools.img \
		vcf-concat {input} > {output}
		"""

def concat_vcf_stats(wildcards):
	return expand("/scratch/n/nicholsa/zyfniu/AN_WGS/results/mutect2/" + wildcards.tumor + "_vs_" + wildcards.patient + "-N" + "/unfiltered_" + wildcards.tumor + "_vs_" + wildcards.patient + "-N.{chr}.vcf.stats", chr = CHROMOSOMES)

rule merge_stats:
	input:
		concat_vcf_stats
	output:
		"/scratch/n/nicholsa/zyfniu/AN_WGS/results/mutect2/{tumor}_vs_{patient}-N/unfiltered_{tumor}_vs_{patient}-N_merged.stats"
	threads: 2
	group: "variantCalling"
	run:
		import os, glob
		cmd = "singularity exec -B $SCRATCH/igenomes_ref,$SCRATCH/HPV_WGS/raw /gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/gatk-4.1.8.img gatk MergeMutectStats "
		for input_file in input:
			cmd = cmd + " --stats " + input_file
		cmd = cmd + " -O {output}"
		shell(cmd)

def concat_vcf_f12(wildcards):
	return expand("/scratch/n/nicholsa/zyfniu/AN_WGS/results/mutect2/" + wildcards.tumor + "_vs_" + wildcards.patient + "-N" + "/unfiltered_" + wildcards.tumor + "_vs_" + wildcards.patient + "-N_f12.{chr}.tar.gz", chr = CHROMOSOMES)

rule gatk_LearnOrientationModel:
	input:
		concat_vcf_f12
	output:
		model = "/scratch/n/nicholsa/zyfniu/AN_WGS/results/mutect2/{tumor}_vs_{patient}-N/unfiltered_{tumor}_read_orientation_model.tar.gz"
	group: "variantCalling"
	run:
		import os, glob
		cmd = "singularity exec -B $SCRATCH/igenomes_ref,$SCRATCH/HPV_WGS/raw /gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/gatk-4.1.8.img gatk LearnReadOrientationModel "
		for input_file in input:
			cmd = cmd + " -I " + input_file
		cmd = cmd + " -O {output}"
		shell(cmd)

rule gatk_get_pileupsummaries_normal:
	input:
		##Downloaded from https://console.cloud.google.com/storage/browser/_details/gatk-best-practices/somatic-hg38/
		common_biallelic_vcf =  REF_exac_common,
		normal = "/scratch/n/nicholsa/zyfniu/AN_WGS/orphan/{patient}-N/Recal/{patient}-N.recal.bam"
	output:
		summary = "/scratch/n/nicholsa/zyfniu/AN_WGS/results/mutect2/{tumor}_vs_{patient}-N/{patient}-N_pileupsummaries.table"
	threads: 2
	group: "variantCalling"
	shell:
		"""
		singularity exec -B $SCRATCH/igenomes_ref,$SCRATCH/HPV_WGS/raw /gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/gatk-4.1.8.img \
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
		primary = "/scratch/n/nicholsa/zyfniu/AN_WGS/orphan/{tumor}/Recal/{tumor}.recal.bam"
	output:
		summary = "/scratch/n/nicholsa/zyfniu/AN_WGS/results/mutect2/{tumor}_vs_{patient}-N/{tumor}_pileupsummaries.table"
	threads: 2
	group: "variantCalling"
	shell:
		"""
		singularity exec -B $SCRATCH/igenomes_ref,$SCRATCH/HPV_WGS/raw /gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/gatk-4.1.8.img \
		gatk --java-options "-Xmx8g" GetPileupSummaries \
		-I {input.primary} \
		-V {input.common_biallelic_vcf} \
		-L {input.common_biallelic_vcf} \
		-O {output.summary}
		"""

rule gatk_calcContam_primary:
	input:
		normal = "/scratch/n/nicholsa/zyfniu/AN_WGS/results/mutect2/{tumor}_vs_{patient}-N/{patient}-N_pileupsummaries.table",
		primary = "/scratch/n/nicholsa/zyfniu/AN_WGS/results/mutect2/{tumor}_vs_{patient}-N/{tumor}_pileupsummaries.table"
	output:
		contamTable = temp("/scratch/n/nicholsa/zyfniu/AN_WGS/results/mutect2/{tumor}_vs_{patient}-N/{tumor}_calContam.table"),
		segment = temp("/scratch/n/nicholsa/zyfniu/AN_WGS/results/mutect2/{tumor}_vs_{patient}-N/{tumor}_tumor.segment")
	threads: 2
	group: "variantCalling"
	shell:
		"""
		singularity exec -B $SCRATCH/igenomes_ref,$SCRATCH/HPV_WGS/raw /gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/gatk-4.1.8.img \
		gatk --java-options "-Xmx8g" CalculateContamination \
		-I {input.normal} \
		-matched {input.primary} \
		--tumor-segmentation {output.segment} \
		-O {output.contamTable}
		"""

rule index_unfiltered_vcf:
	input:
		vcf = "/scratch/n/nicholsa/zyfniu/AN_WGS/results/mutect2/{tumor}_vs_{patient}-N/unfiltered_{tumor}_vs_{patient}-N.vcf",
	output:
		vcf_tbi = "/scratch/n/nicholsa/zyfniu/AN_WGS/results/mutect2/{tumor}_vs_{patient}-N/unfiltered_{tumor}_vs_{patient}-N.vcf.idx"
	threads: 2
	group: "variantCalling"
	shell:
		"""
		singularity exec -B $SCRATCH/igenomes_ref,$SCRATCH/HPV_WGS/raw /gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/gatk-4.1.8.img \
		gatk --java-options "-Xmx8g" IndexFeatureFile \
		-I {input.vcf} \
		--output {output.vcf_tbi}
		"""

rule gatk_filterMutect:
	input:
		vcf_tbi = "/scratch/n/nicholsa/zyfniu/AN_WGS/results/mutect2/{tumor}_vs_{patient}-N/unfiltered_{tumor}_vs_{patient}-N.vcf.idx",
		vcf = "/scratch/n/nicholsa/zyfniu/AN_WGS/results/mutect2/{tumor}_vs_{patient}-N/unfiltered_{tumor}_vs_{patient}-N.vcf",
		model = "/scratch/n/nicholsa/zyfniu/AN_WGS/results/mutect2/{tumor}_vs_{patient}-N/unfiltered_{tumor}_read_orientation_model.tar.gz",
		stats = "/scratch/n/nicholsa/zyfniu/AN_WGS/results/mutect2/{tumor}_vs_{patient}-N/unfiltered_{tumor}_vs_{patient}-N_merged.stats"
	output:
		filter_vcf = temp("/scratch/n/nicholsa/zyfniu/AN_WGS/results/mutect2/{tumor}_vs_{patient}-N/filtered_{tumor}_vs_{patient}-N.{chr}.vcf")
	threads: 2
	group: "variantCalling"
	shell:
		"""
		singularity exec -B $SCRATCH/igenomes_ref,$SCRATCH/HPV_WGS/raw /gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/gatk-4.1.8.img \
		gatk --java-options "-Xmx8g" FilterMutectCalls \
		--ob-priors {input.model} \
		-stats {input.stats} \
		-R {REF_fasta} \
		-V {input.vcf} \
		-L {wildcards.chr} \
		--output {output.filter_vcf}
		"""

def concat_vcf_filtered(wildcards):
	return	expand("/scratch/n/nicholsa/zyfniu/AN_WGS/results/mutect2/" + wildcards.tumor + "_vs_" + wildcards.patient + "-N" + "/filtered_" + wildcards.tumor + "_vs_" + wildcards.patient + "-N.{chr}.vcf", chr = CHROMOSOMES)

rule merge_mutect2_vcf_filtered:
	input:
		concat_vcf_filtered
	output:
		"/scratch/n/nicholsa/zyfniu/AN_WGS/results/mutect2/{tumor}_vs_{patient}-N/filtered_{tumor}_vs_{patient}-N.vcf"
	threads: 2
	group: "variantCalling"
	shell:
		"singularity exec -B $SCRATCH/igenomes_ref,$SCRATCH/HPV_WGS/raw /gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/vcftools.img vcf-concat {input} > {output}"


rule index_filtered_vcf:
	input:
		vcf = "/scratch/n/nicholsa/zyfniu/AN_WGS/results/mutect2/{tumor}_vs_{patient}-N/filtered_{tumor}_vs_{patient}-N.vcf",
	output:
		vcf_tbi = "/scratch/n/nicholsa/zyfniu/AN_WGS/results/mutect2/{tumor}_vs_{patient}-N/filtered_{tumor}_vs_{patient}-N.vcf.idx"
	threads: 2
	group: "variantCalling"
	shell:
		"""
		singularity exec -B $SCRATCH/igenomes_ref,$SCRATCH/HPV_WGS/raw /gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/gatk-4.1.8.img \
		gatk --java-options "-Xmx8g" IndexFeatureFile \
		-I {input.vcf} \
		--output {output.vcf_tbi}
		"""

rule config_manta:
	input:
			normal = "/scratch/n/nicholsa/zyfniu/AN_WGS/orphan/{patient}-N/Recal/{patient}-N.recal.bam",
			tumor = "/scratch/n/nicholsa/zyfniu/AN_WGS/orphan/{tumor}/Recal/{tumor}.recal.bam"
	output: temp("/scratch/n/nicholsa/zyfniu/AN_WGS/temp/Manta/{tumor}_vs_{patient}-N/runWorkflow.py")
	threads: 2
	group: "variantCalling"
	shell:
		"""
		singularity exec -B $SCRATCH/igenomes_ref,$SCRATCH/HPV_WGS/raw /gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/nfcore-sarek-2.6.img \
		configManta.py \
		--normalBam {input.normal} \
		--tumorBam  {input.tumor} \
		--reference {REF_fasta} \
		--runDir temp/Manta/{wildcards.tumor}_vs_{wildcards.patient}-N/{wildcards.tumor}
		"""

rule manta:
	input:
			normal = "/scratch/n/nicholsa/zyfniu/AN_WGS/orphan/{patient}-N/Recal/{patient}-N.recal.bam",
			tumor = "/scratch/n/nicholsa/zyfniu/AN_WGS/orphan/{tumor}/Recal/{tumor}.recal.bam",
			script = "/scratch/n/nicholsa/zyfniu/AN_WGS/temp/Manta/{tumor}_vs_{patient}-N/runWorkflow.py"
	output:
			sv = temp("/scratch/n/nicholsa/zyfniu/AN_WGS/results/Manta/{tumor}_vs_{patient}-N/Manta_{tumor}_vs_{patient}-N.candidateSV.vcf.gz"),
			smallindel = temp("/scratch/n/nicholsa/zyfniu/AN_WGS/results/Manta/{tumor}_vs_{patient}-N/Manta_{tumor}_vs_{patient}-N.candidateSmallIndels.vcf.gz"),
			diploidSV = temp("/scratch/n/nicholsa/zyfniu/AN_WGS/results/Manta/{tumor}_vs_{patient}-N/Manta_{tumor}_vs_{patient}-N.diploidSV.vcf.gz"),
			somaticSV = temp("/scratch/n/nicholsa/zyfniu/AN_WGS/results/Manta/{tumor}_vs_{patient}-N/Manta_{tumor}_vs_{patient}-N.somaticSV.vcf.gz")
	threads: 80
	group: "variantCalling"
	shell:
		"""
		singularity exec -B $SCRATCH/igenomes_ref,$SCRATCH/HPV_WGS/raw /gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/nfcore-sarek-2.6.img \
		python {input.script} -m local -j {threads}
		mv temp/Manta/{wildcards.tumor}_vs_{wildcards.patient}-N/results/variants/candidateSmallIndels.vcf.gz \
		results/Manta/{wildcards.tumor}_vs_{wildcards.patient}-N/Manta_{wildcards.tumor}_vs_{wildcards.patient}-N.candidateSmallIndels.vcf.gz
		mv temp/Manta/{wildcards.tumor}_vs_{wildcards.patient}-N/results/variants/candidateSmallIndels.vcf.gz.tbi \
		results/Manta/{wildcards.tumor}_vs_{wildcards.patient}-N/Manta_{wildcards.tumor}_vs_{wildcards.patient}-N.candidateSmallIndels.vcf.gz.tbi
		mv temp/Manta/{wildcards.tumor}_vs_{wildcards.patient}-N/results/variants/candidateSV.vcf.gz \
		results/Manta/{wildcards.tumor}_vs_{wildcards.patient}-N/Manta_{wildcards.tumor}_vs_{wildcards.patient}-N.candidateSV.vcf.gz
		mv temp/Manta/{wildcards.tumor}_vs_{wildcards.patient}-N/results/variants/candidateSV.vcf.gz.tbi \
		results/Manta/{wildcards.tumor}_vs_{wildcards.patient}-N/Manta_{wildcards.tumor}_vs_{wildcards.patient}-N.candidateSV.vcf.gz.tbi
		mv temp/Manta/{wildcards.tumor}_vs_{wildcards.patient}-N/results/variants/diploidSV.vcf.gz \
		results/Manta/{wildcards.tumor}_vs_{wildcards.patient}-N/Manta_{wildcards.tumor}_vs_{wildcards.patient}-N.diploidSV.vcf.gz
		mv temp/Manta/{wildcards.tumor}_vs_{wildcards.patient}-N/results/variants/diploidSV.vcf.gz.tbi \
		results/Manta/{wildcards.tumor}_vs_{wildcards.patient}-N/Manta_{wildcards.tumor}_vs_{wildcards.patient}-N.diploidSV.vcf.gz.tbi
		mv temp/Manta/{wildcards.tumor}_vs_{wildcards.patient}-N/results/variants/somaticSV.vcf.gz \
		results/Manta/{wildcards.tumor}_vs_{wildcards.patient}-N/Manta_{wildcards.tumor}_vs_{wildcards.patient}-N.somaticSV.vcf.gz
		mv temp/Manta/{wildcards.tumor}_vs_{wildcards.patient}-N/results/variants/somaticSV.vcf.gz.tbi \
		results/Manta/{wildcards.tumor}_vs_{wildcards.patient}-N/Manta_{wildcards.tumor}_vs_{wildcards.patient}-N.somaticSV.vcf.gz.tbi
		"""

rule annotate_mutect2:
	input:
		vcf = "/scratch/n/nicholsa/zyfniu/AN_WGS/results/mutect2/{tumor}_vs_{patient}-N/filtered_{tumor}_vs_{patient}-N.vcf"
	output:
			annotatedvcf = temp("/scratch/n/nicholsa/zyfniu/AN_WGS/results/mutect2/{tumor}_vs_{patient}-N/{tumor}_vs_{patient}-N_snpEff.ann.vcf")
	threads: 2
	group: "variantCalling"
	shell:
		"""
		cd results/mutect2/{wildcards.tumor}
		singularity exec -B $SCRATCH/AN_WGS $SCRATCH/singularity_images/nfcore-sareksnpeff-2.6.GRCh38.img \
		snpEff -Xmx8g \
		GRCh38.86 \
		-csvStats {wildcards.tumor}_snpEff.csv \
		-nodownload \
		-canon \
		{input} \
		> {output}
		mv snpEff_summary.html {wildcards.tumor}_snpEff.html
		"""

rule zip_snpeff:
	input:
		annotatedvcf = "/scratch/n/nicholsa/zyfniu/AN_WGS/results/mutect2/{tumor}_vs_{patient}-N/{tumor}_vs_{patient}-N_snpEff.ann.vcf"
	output: annotatedvcf = "/scratch/n/nicholsa/zyfniu/AN_WGS/results/mutect2/{tumor}_vs_{patient}-N/{tumor}_vs_{patient}-N_snpEff.ann.vcf.gz"
	threads: 2
	group: "variantCalling"
	shell:
		"""
		singularity exec /gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/nfcore-sarek-2.6.img bgzip < {input} > {output}
		singularity exec /gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/nfcore-sarek-2.6.img tabix {output}
		"""

rule annotate_manta:
	input:
		vcf = "/scratch/n/nicholsa/zyfniu/AN_WGS/results/Manta/{tumor}_vs_{patient}-N/Manta_{tumor}_vs_{patient}-N.{structure}.vcf.gz"
	output:
		annotatedvcf = temp("/scratch/n/nicholsa/zyfniu/AN_WGS/results/Manta/{tumor}_vs_{patient}-N/Manta_snpeff_{tumor}_vs_{patient}-N.{structure}.ann.vcf")
	threads: 80
	group: "variantCalling"
	shell:
		"""
		cd results/mutect2/{wildcards.tumor}/
		singularity exec -B $SCRATCH/AN_WGS $SCRATCH/singularity_images/nfcore-sareksnpeff-2.6.GRCh38.img \
		snpEff -Xmx8g \
		GRCh38.86 \
		-csvStats {wildcards.tumor}_{wildcards.structure}_snpEff.csv \
		-nodownload \
		-canon \
		{input} \
		> {output}
		mv snpEff_summary.html {wildcards.tumor}_{wildcards.structure}_snpEff.html
		"""

rule zip_manta:
	input:
		annotatedvcf = "/scratch/n/nicholsa/zyfniu/AN_WGS/results/Manta/{tumor}_vs_{patient}-N/Manta_snpeff_{tumor}_vs_{patient}-N.{structure}.ann.vcf"
	output: annotatedvcf = "/scratch/n/nicholsa/zyfniu/AN_WGS/results/Manta/{tumor}_vs_{patient}-N/Manta_snpeff_{tumor}_vs_{patient}-N.{structure}.ann.vcf.gz"
	threads: 2
	group: "variantCalling"
	shell:
		"""
		singularity exec /gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/nfcore-sarek-2.6.img bgzip < {input} > {output}
		singularity exec /gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/nfcore-sarek-2.6.img tabix {output}
		"""

###
### COPY NUMBER USING ASCAT
###
rule alleleCount:
	input:
		bam = "orphan/{sample}/Recal/{sample}.recal.bam",
		index = "orphan/{sample}/Recal/{sample}.recal.bai",
		acloci = "/scratch/n/nicholsa/zyfniu/igenomes_ref/1000G_phase3_GRCh38_maf0.3.loci"
	output: "/scratch/n/nicholsa/zyfniu/AN_WGS/results/ASCAT/alleleCount/{sample}.alleleCount"
	threads: 2
	group: "variantCalling"
	shell:
		"""
		singularity exec  -B $SCRATCH/igenomes_ref /gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/nfcore-sarek-2.6.img \
		alleleCounter \
		-l {input.acloci} \
		-r {REF_fasta} \
		-b {input.bam} \
		-o results/ASCAT/alleleCount/{wildcards.sample}.alleleCount
		"""

def getGender(wildcards):
		return expand("{gender}",gender = INPUT[INPUT.Patient == wildcards.patient].Sex.drop_duplicates())

###convert allele script from https://bitbucket.org/malinlarsson/somatic_wgs_pipeline/src/master/convertAlleleCounts.r
rule ConvertAlleleCounts:
	input:
		normal = "/scratch/n/nicholsa/zyfniu/AN_WGS/results/ASCAT/alleleCount/{patient}-N.alleleCount",
		tumor = "/scratch/n/nicholsa/zyfniu/AN_WGS/results/ASCAT/alleleCount/{tumor}.alleleCount",
		script = "/scratch/n/nicholsa/zyfniu/AN_WGS/AN_WGS/convertAlleleCounts.r"
	output:
		normalBaf = "/scratch/n/nicholsa/zyfniu/AN_WGS/results/ASCAT/{tumor}_vs_{patient}-N/{patient}-N.BAF",
		tumorBaf = "/scratch/n/nicholsa/zyfniu/AN_WGS/results/ASCAT/{tumor}_vs_{patient}-N/{tumor}.BAF",
		normalLogr = "/scratch/n/nicholsa/zyfniu/AN_WGS/results/ASCAT/{tumor}_vs_{patient}-N/{patient}-N.LogR",
		tumorLogr = "/scratch/n/nicholsa/zyfniu/AN_WGS/results/ASCAT/{tumor}_vs_{patient}-N/{tumor}.LogR"
	threads: 2
	params: gender = getGender
	group: "variantCalling"
	shell:
		"""
		cd results/ASCAT/{wildcards.tumor}_vs_{wildcards.patient}-N
		singularity exec -B $SCRATCH/AN_WGS /gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/nfcore-sarek-2.6.img Rscript \
		/scratch/n/nicholsa/zyfniu/AN_WGS/AN_WGS/convertAlleleCounts.r \
		{wildcards.tumor} /scratch/n/nicholsa/zyfniu/AN_WGS/results/alleleCount/{wildcards.tumor}.alleleCount \
		{wildcards.patient}-N /scratch/n/nicholsa/zyfniu/AN_WGS/results/alleleCount/{wildcards.patient}-N.alleleCount \
		{params.gender}
		"""

rule ascat:
	input:
		normalBaf = "/scratch/n/nicholsa/zyfniu/AN_WGS/results/ASCAT/{tumor}_vs_{patient}-N/{patient}-N.BAF",
		tumorBaf = "/scratch/n/nicholsa/zyfniu/AN_WGS/results/ASCAT/{tumor}_vs_{patient}-N/{tumor}.BAF",
		normalLogr = "/scratch/n/nicholsa/zyfniu/AN_WGS/results/ASCAT/{tumor}_vs_{patient}-N/{patient}-N.LogR",
		tumorLogr = "/scratch/n/nicholsa/zyfniu/AN_WGS/results/ASCAT/{tumor}_vs_{patient}-N/{tumor}.LogR",
		acLociGC = "/scratch/n/nicholsa/zyfniu/igenomes_ref/1000G_phase3_GRCh38_maf0.3.loci.gc"
	output:
		results = "/scratch/n/nicholsa/zyfniu/AN_WGS/results/ASCAT/{tumor}_vs_{patient}-N/{tumor}_vs_{patient}-N.tumor.cnvs.txt"
	params:
		gender = getGender
	threads: 8
	group: "variantCalling"
	shell:
		"""
		for f in /scratch/n/nicholsa/zyfniu/AN_WGS/results/ASCAT/{wildcards.tumor}_vs_{wildcards.patient}-N/*BAF \
		/scratch/n/nicholsa/zyfniu/AN_WGS/results/ASCAT/{wildcards.tumor}_vs_{wildcards.patient}-N/*LogR; \
		do sed \'s/chr//g\' $f > /scratch/n/nicholsa/zyfniu/AN_WGS/results/ASCAT/{wildcards.tumor}_vs_{wildcards.patient}-N/tmpFile; \
		mv /scratch/n/nicholsa/zyfniu/AN_WGS/results/ASCAT/{wildcards.tumor}_vs_{wildcards.patient}-N/tmpFile $f;done
		cd results/ASCAT/{wildcards.tumor}_vs_{wildcards.patient}-N
		singularity exec -B $SCRATCH/AN_WGS,$SCRATCH/igenomes_ref /gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/nfcore-sarek-2.6.img \
		Rscript /scratch/n/nicholsa/zyfniu/AN_WGS/AN_WGS/run_ascat.r \
        --tumorbaf {input.tumorBaf} \
        --tumorlogr {input.tumorLogr} \
        --normalbaf {input.normalBaf} \
        --normallogr {input.normalLogr} \
        --tumorname {wildcards.tumor} \
        --basedir results/ASCAT/{wildcards.tumor}_vs_{wildcards.patient}-N \
        --gcfile {input.acLociGC} \
        --gender {params.gender}
		cd ../../../
		mv results/ASCAT/{wildcards.tumor}_vs_{wildcards.patient}-N/{wildcards.tumor}.cnvs.txt {output.results}
		"""
###
### RUN MULTIQC
###

# def getFinishedMutect2(wildcards):
# 		return expand("results/mutect2/{tumor}/{tumor}_vs_{patient}-N_snpEff.ann.vcf.gz",
# 		tumor = INPUT[(INPUT.n_vs_t == 1)&(INPUT.Patient == wildcards.patient)].Sample.drop_duplicates(),
# 		patient = INPUT[(INPUT.n_vs_t == 0)&(INPUT.Patient == wildcards.patient)].Patient.drop_duplicates())
#
# def getFinishedManta(wildcards):
# 		return expand("results/Manta/{tumor}/Manta_snpeff_{tumor}_vs_{patient}-N.{structure}.ann.vcf.gz",
# 		tumor = INPUT[(INPUT.n_vs_t == 1)&(INPUT.Patient == wildcards.patient)].Sample.drop_duplicates(),
# 		patient = INPUT[(INPUT.n_vs_t == 0)&(INPUT.Patient == wildcards.patient)].Patient.drop_duplicates(),
# 		structure = ["candidateSV","candidateSmallIndels","diploidSV","somaticSV"])
#
# def getFinishedASCAT(wildcards):
# 		return expand("results/ASCAT/{tumor}_vs_{patient}-N/{tumor}_vs_{patient}-N.tumor.cnvs.txt",
# 		tumor = INPUT[(INPUT.n_vs_t == 1)&(INPUT.Patient == wildcards.patient)].Sample.drop_duplicates(),
# 		patient = INPUT[(INPUT.n_vs_t == 0)&(INPUT.Patient == wildcards.patient)].Patient.drop_duplicates())

# rule multiqc:
# 	input:
# 		getFinishedMutect2,
# 		getFinishedManta,
# 		getFinishedASCAT
# 	output:
# 		"reports/{patient}_multiqc.html"
# 	group: "multiqc"
# 	threads: 80
# 	shell:
# 		"""
# 		singularity exec $SCRATCH/singularity_images/nfcore-sarek-2.6.img multiqc . -n {output} \
# 		results/mutect2/{wildcards.patient} results/ASCAT/*{wildcards.patient}* results/Manta/{wildcards.patient}
# 		"""

###
###	END
###
