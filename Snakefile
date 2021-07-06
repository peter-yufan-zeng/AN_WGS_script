#### USAGE
#### 1) First load snakemake
#### module load NiaEnv/2018a
#### module load python/3.6.4-anaconda5.1.0
#### source activate snakemake
#### snakemake -s AN_WGS_script/Snakefile  --cores 1 -j 40 --cluster "sbatch -N 1 -t 20:00:00 --ntasks 80 --output=logs/%x_%j.log" --ri \
#### --config input=AN_WGS_script/Sample/20200907.tsv.txt outdir=$SCRATCH/AN_WGS_script/20200908_HPV_HNSCC_WGS

import pandas as pd
import os

#include: "snakemake_scripts/hla.snk"
#localrules: all

###
### LOAD SAMPLES
###
SCRATCH = "/gpfs/fs0/scratch/n/nicholsa/zyfniu"
print("\n***INPUT FILE: " + config['input'] + "***\n")
INPUT = pd.read_table(config['input'],names = ['Patient','Sex','n_vs_t','Sample','Lane','Fastq1','Fastq2'])
INPUT['Lane'] = INPUT.Lane.apply(str)
INPUT['Sample_Lane'] = INPUT.Sample + "_" + INPUT.Lane
SAMPLE =  INPUT['Sample'].unique()
TMP = "temp"
PAT = INPUT.Patient.drop_duplicates()
TUMOR = INPUT[(INPUT.n_vs_t == 1)].Sample.drop_duplicates()
SAMPLE_LANE = INPUT.Sample_Lane
print(INPUT)


### PRINT OUTPUT DIRECTORY
OUTDIR = config['outdir']
print("***OUTPUT DIRECTORY: " + OUTDIR + "***")

###
### GET SAMPLES TO USE MULTI-TUMOUR VARIANT CALLING
###
Multi = INPUT[(INPUT.n_vs_t == 1)]
Multi = Multi[['Patient','Sample']].drop_duplicates()
Multi.loc[:,'Combined'] = ""
Multi.loc[:,"num_tumours"] = 1
for x in Multi['Patient'].drop_duplicates():
	COMBINED_mutect_samples = ""
	print(x)
	n = 1
	for y in Multi[Multi.Patient == x].Sample.drop_duplicates():
		print(y)
		COMBINED_mutect_samples = y + "_" + COMBINED_mutect_samples
		Multi.loc[(Multi.Patient == x),'num_tumours'] = n
		n = n + 1
	COMBINED_mutect_samples = COMBINED_mutect_samples[:-1]
	Multi.loc[(Multi.Patient == x),'Combined']= COMBINED_mutect_samples

Multi['path'] = OUTDIR + "/results/mutect2/all_"  + Multi.Combined + "_vs_" + Multi.Patient + "-N/" \
+ Multi.Combined + "_vs_" + Multi.Patient + "-N_snpEff.ann.vcf.gz"
Multi = Multi[Multi.num_tumours > 1]
#Multi = Multi[['Patient','Combined',"path"]].drop_duplicates()



###
###	REFERENCE FILES
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


###
### Final Results
###
rule all:
	input:
		###FASTQC + BWA Alignment
		expand(OUTDIR + "/QC/{sample_lane}/{sample_lane}_R1_fastqc.html", sample_lane = SAMPLE_LANE),
		expand(OUTDIR + "/Recal/{sample}.recal.bam", sample = SAMPLE),
		### using zip https://endrebak.gitbooks.io/the-snakemake-book/chapters/expand/expand.html
		###BAMQC + samtools_stats
		expand(OUTDIR + "/QC/{sample}/{sample}.samtools.stats.out",sample = SAMPLE),
		expand(OUTDIR + "/QC/{sample}/bamQC/qualimapReport.html",sample = SAMPLE),
		####
		####VARIANT CALLING OUTPUTS
		####
		expand(OUTDIR + "/results/mutect2/{tumor}_vs_{patient}-N/{tumor}_vs_{patient}-N_snpEff.ann.vcf.gz",zip, patient = [x.rsplit('-',1)[0] for x in TUMOR],tumor = TUMOR),
		expand(OUTDIR + "/results/mutect2/{tumor}_vs_{patient}-N/{tumor}_vs_{patient}-N_snpEff.ann.passOnly.vcf.gz",zip, patient = [x.rsplit('-',1)[0] for x in TUMOR],tumor = TUMOR),
		Multi.path.drop_duplicates().tolist(),
		#### Manta
		expand(OUTDIR + "/results/Manta/{tumor}_vs_{patient}-N/Manta_snpeff_{tumor}_vs_{patient}-N.candidateSV.ann.vcf.gz",zip, patient = [x.rsplit('-',1)[0] for x in TUMOR],tumor = TUMOR),
		expand(OUTDIR + "/results/Manta/{tumor}_vs_{patient}-N/Manta_snpeff_{tumor}_vs_{patient}-N.candidateSmallIndels.ann.vcf.gz",zip, patient = [x.rsplit('-',1)[0] for x in TUMOR],tumor = TUMOR),
		expand(OUTDIR + "/results/Manta/{tumor}_vs_{patient}-N/Manta_snpeff_{tumor}_vs_{patient}-N.diploidSV.ann.vcf.gz",zip, patient = [x.rsplit('-',1)[0] for x in TUMOR],tumor = TUMOR),
		expand(OUTDIR + "/results/Manta/{tumor}_vs_{patient}-N/Manta_snpeff_{tumor}_vs_{patient}-N.somaticSV.ann.vcf.gz",zip, patient = [x.rsplit('-',1)[0] for x in TUMOR],tumor = TUMOR),
		#### ASCAT
		expand(OUTDIR + "/results/ASCAT/{tumor}_vs_{patient}-N/{tumor}_vs_{patient}-N.tumor.cnvs.txt",zip, tumor = TUMOR, patient = [x.rsplit('-',1)[0] for x in TUMOR])
		#### QC
		#expand("reports/{patient}_multiqc.html",patient = PAT)
	threads: 80
	run:
		import time
		os.system("singularity exec " + OUTDIR + " $SCRATCH/singularity_images/nfcore-sarek-2.6.img multiqc " + OUTDIR + " -n  " + OUTDIR + "/reports/multiqc_" + time.strftime('%Y%m%d') + ".html")
		os.system("mv multiqc_data* " + OUTDIR + "/QC/")
		os.system("cp " + config['input'] + " " + OUTDIR +  "/Sample_" + time.strftime('%Y%m%d') + ".tsv")
		os.system("rm " + OUTDIR + "/temp -rf")
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
		o1 = OUTDIR + "/QC/{sample_lane}/{sample_lane}_R1_fastqc.html",
		o2 = OUTDIR + "/QC/{sample_lane}/{sample_lane}_R2_fastqc.html",
		o3 = OUTDIR + "/QC/{sample_lane}/{sample_lane}_R1_fastqc.zip",
		o4 = OUTDIR + "/QC/{sample_lane}/{sample_lane}_R2_fastqc.zip"
	threads: 80
	group: "align"
	shell:
		"""
		singularity exec -B $SCRATCH/igenomes_ref,$SCRATCH/AN_WGS/raw,{OUTDIR} /gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/nfcore-sarek-2.6.img \
		fastqc -t {threads} {input} --outdir={OUTDIR}/QC/{wildcards.sample_lane}/
		mv {OUTDIR}/QC/{wildcards.sample_lane}/*R1*fastqc.html {OUTDIR}/QC/{wildcards.sample_lane}/{wildcards.sample_lane}_R1_fastqc.html
		mv {OUTDIR}/QC/{wildcards.sample_lane}/*R2*fastqc.html {OUTDIR}/QC/{wildcards.sample_lane}/{wildcards.sample_lane}_R2_fastqc.html
		mv {OUTDIR}/QC/{wildcards.sample_lane}/*R1*fastqc.zip {OUTDIR}/QC/{wildcards.sample_lane}/{wildcards.sample_lane}_R1_fastqc.zip
		mv {OUTDIR}/QC/{wildcards.sample_lane}/*R2*fastqc.zip {OUTDIR}/QC/{wildcards.sample_lane}/{wildcards.sample_lane}_R2_fastqc.zip
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
		temp(OUTDIR +"/orphan/{sample_lane}/{sample_lane}.bwa.sam")
	threads: 80
	group: "align"
	params:
		readGroup = createRG
	shell:
		"""
		singularity exec -B $SCRATCH/igenomes_ref,$SCRATCH/AN_WGS/raw /gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/nfcore-sarek-2.6.img \
		bwa mem -K 100000000 -R \"{params.readGroup}\" -B 3 -t {threads} -M {REF_fasta} \
			{input.fastq1} {input.fastq2} > {output}
		"""

rule sort_sam_to_bam:
	input:
		sam = OUTDIR +"/orphan/{sample_lane}/{sample_lane}.bwa.sam"
	output:
		bam = OUTDIR +"/orphan/{sample_lane}/{sample_lane}.bwa.bam"
	threads: 70
	params:
		temp = OUTDIR +"/orphan/{sample_lane}/"
	group: "align"
	shell:
		"""
		singularity exec -B $SCRATCH/igenomes_ref,$SCRATCH/AN_WGS/raw,{OUTDIR} /gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/nfcore-sarek-2.6.img \
		samtools sort -T {params.temp} --threads {threads} -m 2G {input.sam} > {output}
		"""

def get_bams_to_merge(wildcards):
	SAMPLE_LANE = INPUT[INPUT.Sample == wildcards.sample].Sample_Lane
	return expand(OUTDIR +"/orphan/{sample_lane}/{sample_lane}.bwa.bam", sample_lane = SAMPLE_LANE)

rule merge_bam_mapped_and_index:
	input:
		get_bams_to_merge
	output:
		temp(OUTDIR +"/orphan/{sample}/{sample}.merged.bam")
	threads: 10
	group: "merge_markduplicate"
	run:
		if len(input) == 1:
			shell("mv {input} {output}")
			shell("singularity exec -B $SCRATCH/igenomes_ref,$SCRATCH/AN_WGS/raw,{OUTDIR} /gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/nfcore-sarek-2.6.img samtools index -@ {threads} {output}")
		elif len(input) > 1:
			#shell("mv {input} {output}")
			shell("singularity exec -B $SCRATCH/igenomes_ref,$SCRATCH/AN_WGS/raw,{OUTDIR} /gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/nfcore-sarek-2.6.img samtools merge --threads {threads} - {input} | tee {output} | singularity exec -B $SCRATCH/igenomes_ref,$SCRATCH/AN_WGS/raw,{OUTDIR} /gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/nfcore-sarek-2.6.img samtools index -@ {threads} - ")

rule markdup:
	input:
		OUTDIR +"/orphan/{sample}/{sample}.merged.bam"
	output:
		bam = temp(OUTDIR +"/orphan/MD/{sample}.md.bam"),
		metric = OUTDIR +"/orphan/MD/{sample}.bam.metric"
	threads: 20
	group: "merge_markduplicate"
	shell:
		"""
		singularity exec -B $SCRATCH/igenomes_ref,$SCRATCH/AN_WGS,$SCRATCH/AN_WGS/raw,{OUTDIR} /gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/gatk-4.1.8.img \
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
		OUTDIR +"/orphan/MD/{sample}.md.bam"
	output:
		temp(OUTDIR +"/orphan/{sample}/Recal/{sample}.{chr}.recal.table")
	threads: 2
	group: "recalibrator"
	shell:
		"""
		singularity exec -B $SCRATCH/igenomes_ref,$SCRATCH/AN_WGS,$SCRATCH/AN_WGS/raw,{OUTDIR} /gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/gatk-4.1.8.img \
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
	return expand(OUTDIR +"/orphan/" + wildcards.sample + "/Recal/" + wildcards.sample + ".{chr}.recal.table",chr = CHROMOSOMES)

rule gather_recal_tbl:
	input:
		recal_tbl_to_gather
	output:
		temp(OUTDIR +"/orphan/{sample}/Recal/{sample}.recal.table")
	threads: 2
	group: "recalibrator"
	run:
		command = 'singularity exec -B $SCRATCH/igenomes_ref,$SCRATCH/AN_WGS/raw,{OUTDIR} /gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/nfcore-sarek-2.6.img gatk --java-options \"-Xmx8g\" GatherBQSRReports'
		for i in input:
			command = command + " -I " + i
		command = command + ' -O {output}'
		shell(command)

rule ApplyBQSR:
	input:
		bam = OUTDIR +"/orphan/MD/{sample}.md.bam",
		recal_table = OUTDIR +"/orphan/{sample}/Recal/{sample}.recal.table"
	output:
		temp("orphan/{sample}/Recal/{sample}.recal.{chr}.bam")
	group: "recalibrator"
	shell:
		"""
		singularity exec -B $SCRATCH/igenomes_ref,$SCRATCH/AN_WGS,$SCRATCH/AN_WGS/raw,{OUTDIR} /gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/gatk-4.1.8.img \
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
		bam = OUTDIR +"/Recal/{sample}.recal.bam",
		index = OUTDIR +"/Recal/{sample}.recal.bai"
	group: "recalibrator"
	threads: 20
	shell:
		"""
		singularity exec -B $SCRATCH/igenomes_ref,{OUTDIR} /gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/nfcore-sarek-2.6.img \
		samtools merge --threads {threads} {output.bam} {input}
		singularity exec -B $SCRATCH/igenomes_ref,{OUTDIR} /gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/nfcore-sarek-2.6.img \
		samtools index {output.bam} {output.index}
		"""

###
### QC On BAM FILES
###
rule samtools_stats:
	input:
		bam = OUTDIR +"/Recal/{sample}.recal.bam",
		index = OUTDIR +"/Recal/{sample}.recal.bai"
	output:
		stats = OUTDIR + "/QC/{sample}/{sample}.samtools.stats.out"
	group: "variantCalling"
	threads: 2
	shell:
		"""
		singularity exec -B $SCRATCH/igenomes_ref,{OUTDIR} \
			/gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/nfcore-sarek-2.6.img \
		samtools stats {input.bam} > {output}
		"""

rule bamqc:
	input:
		bam = OUTDIR +"/Recal/{sample}.recal.bam",
		index = OUTDIR +"/Recal/{sample}.recal.bai"
	output:
		stats = OUTDIR + "/QC/{sample}/bamQC/qualimapReport.html"
	group: "variantCalling"
	threads: 40
	shell:
		"""
		singularity exec -B $SCRATCH/igenomes_ref,{OUTDIR} \
			/gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/nfcore-sarek-2.6.img \
		qualimap --java-mem-size=100G \
		bamqc \
		-bam {input.bam} \
		--paint-chromosome-limits \
		--genome-gc-distr HUMAN \
		-nt {threads} \
		-skip-duplicated \
		--skip-dup-mode 0 \
		-outdir {OUTDIR}/QC/{wildcards.sample}/bamQC/ \
		-outformat HTML
		"""

###
### VARIANT CALLING FOR EACH TUMOUR INDIVIDUALLY
###

rule mutect2:
	input:
		normal = OUTDIR +"/Recal/{patient}-N.recal.bam",
		tumor = OUTDIR +"/Recal/{tumor}.recal.bam",
		interval = "/scratch/n/nicholsa/zyfniu/igenomes_ref/interval-files-folder/{num}-scattered.interval_list"
	#	recurrence = "{sample}R/Recal/{sample}R.recal.bam"
	output:
		vcf = temp(OUTDIR + "/results/mutect2/{tumor}_vs_{patient}-N/unfiltered_{tumor}_vs_{patient}-N.{num}.vcf"),
		stats = temp(OUTDIR + "/results/mutect2/{tumor}_vs_{patient}-N/unfiltered_{tumor}_vs_{patient}-N.{num}.vcf.stats"),
		f12 = temp(OUTDIR + "/results/mutect2/{tumor}_vs_{patient}-N/unfiltered_{tumor}_vs_{patient}-N_f12.{num}.tar.gz"),
		index =  temp(OUTDIR + "/results/mutect2/{tumor}_vs_{patient}-N/unfiltered_{tumor}_vs_{patient}-N.{num}.vcf.idx")
	threads: 2
	group: "variantCalling"
	shell:
		"""
		singularity exec -B $SCRATCH/igenomes_ref,{OUTDIR} /gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/gatk-4.1.8.img \
		gatk --java-options "-Xmx8g" Mutect2 -R {REF_fasta} \
		-I {input.normal}  \
		-I {input.tumor} \
		-normal {wildcards.patient}-N \
		-L {input.interval} \
		--germline-resource {REF_gnomAD} \
		--panel-of-normals {REF_pon} --f1r2-tar-gz {output.f12} \
		-O {output.vcf}
		"""
##SET NUMBERS FROM 0000 to 01999

NUMBERS = [str(a)+str(b)+str(c)+str(d) for a in range(0,1) for b in range(0,2) for c in range(0,10) for d in range(0,10)]

def concat_vcf(wildcards):
	return expand(OUTDIR + "/results/mutect2/" + wildcards.tumor + "_vs_" + wildcards.patient + "-N" + "/unfiltered_" + wildcards.tumor + "_vs_" + wildcards.patient + "-N.{num}.vcf", num = NUMBERS)

rule merge_mutect2_vcf:
	input:
		concat_vcf
	output:
		temp(OUTDIR + "/results/mutect2/{tumor}_vs_{patient}-N/unfiltered_{tumor}_vs_{patient}-N.vcf")
	threads: 1
	group: "variantCalling"
	shell:
		"""
		singularity exec -B $SCRATCH/igenomes_ref,{OUTDIR} /gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/vcftools.img \
		vcf-concat {input} > {output}
		"""

def concat_vcf_stats(wildcards):
	return expand(OUTDIR + "/results/mutect2/" + wildcards.tumor + "_vs_" + wildcards.patient + "-N" + "/unfiltered_" + wildcards.tumor + "_vs_" + wildcards.patient + "-N.{num}.vcf.stats", num = NUMBERS)

rule merge_stats:
	input:
		concat_vcf_stats
	output:
		temp(OUTDIR + "/results/mutect2/{tumor}_vs_{patient}-N/unfiltered_{tumor}_vs_{patient}-N_merged.stats")
	threads: 2
	group: "variantCalling"
	run:
		import os, glob
		cmd = "singularity exec -B $SCRATCH/igenomes_ref,{OUTDIR} /gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/gatk-4.1.8.img gatk MergeMutectStats "
		for input_file in input:
			cmd = cmd + " --stats " + input_file
		cmd = cmd + " -O {output}"
		shell(cmd)

def concat_vcf_f12(wildcards):
	return expand(OUTDIR + "/results/mutect2/" + wildcards.tumor + "_vs_" + wildcards.patient + "-N" + "/unfiltered_" + wildcards.tumor + "_vs_" + wildcards.patient + "-N_f12.{num}.tar.gz", num = NUMBERS)

rule gatk_LearnOrientationModel:
	input:
		concat_vcf_f12
	output:
		model = OUTDIR + "/results/mutect2/{tumor}_vs_{patient}-N/unfiltered_{tumor}_read_orientation_model.tar.gz"
	group: "variantCalling"
	threads: 2
	run:
		import os, glob
		cmd = "singularity exec -B $SCRATCH/igenomes_ref,{OUTDIR} /gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/gatk-4.1.8.img gatk LearnReadOrientationModel "
		for input_file in input:
			cmd = cmd + " -I " + input_file
		cmd = cmd + " -O {output}"
		shell(cmd)

rule gatk_get_pileupsummaries_normal:
	input:
		##Downloaded from https://console.cloud.google.com/storage/browser/_details/gatk-best-practices/somatic-hg38/
		common_biallelic_vcf =  REF_exac_common,
		normal = OUTDIR + "/orphan/{patient}-N/Recal/{patient}-N.recal.bam"
	output:
		summary = OUTDIR + "/results/mutect2/{tumor}_vs_{patient}-N/{patient}-N_pileupsummaries.table"
	threads: 2
	group: "variantCalling"
	shell:
		"""
		singularity exec -B $SCRATCH/igenomes_ref,{OUTDIR} /gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/gatk-4.1.8.img \
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
		primary = OUTDIR + "/orphan/{tumor}/Recal/{tumor}.recal.bam"
	output:
		summary = OUTDIR + "/results/mutect2/{tumor}_vs_{patient}-N/{tumor}_pileupsummaries.table"
	threads: 2
	group: "variantCalling"
	shell:
		"""
		singularity exec -B $SCRATCH/igenomes_ref,{OUTDIR} /gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/gatk-4.1.8.img \
		gatk --java-options "-Xmx8g" GetPileupSummaries \
		-I {input.primary} \
		-V {input.common_biallelic_vcf} \
		-L {input.common_biallelic_vcf} \
		-O {output.summary}
		"""

rule gatk_calcContam_primary:
	input:
		normal = OUTDIR + "/results/mutect2/{tumor}_vs_{patient}-N/{patient}-N_pileupsummaries.table",
		primary = OUTDIR + "/results/mutect2/{tumor}_vs_{patient}-N/{tumor}_pileupsummaries.table"
	output:
		contamTable = temp(OUTDIR + "/results/mutect2/{tumor}_vs_{patient}-N/{tumor}_calContam.table"),
		segment = temp(OUTDIR + "/results/mutect2/{tumor}_vs_{patient}-N/{tumor}_tumor.segment")
	threads: 2
	group: "variantCalling"
	shell:
		"""
		singularity exec -B $SCRATCH/igenomes_ref,{OUTDIR} /gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/gatk-4.1.8.img \
		gatk --java-options "-Xmx8g" CalculateContamination \
		-I {input.normal} \
		-matched {input.primary} \
		--tumor-segmentation {output.segment} \
		-O {output.contamTable}
		"""

rule index_unfiltered_vcf:
	input:
		vcf = OUTDIR + "/results/mutect2/{tumor}_vs_{patient}-N/unfiltered_{tumor}_vs_{patient}-N.vcf",
	output:
		vcf_tbi = OUTDIR + "/results/mutect2/{tumor}_vs_{patient}-N/unfiltered_{tumor}_vs_{patient}-N.vcf.idx"
	threads: 2
	group: "variantCalling"
	shell:
		"""
		singularity exec -B $SCRATCH/igenomes_ref,{OUTDIR} /gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/gatk-4.1.8.img \
		gatk --java-options "-Xmx8g" IndexFeatureFile \
		-I {input.vcf} \
		--output {output.vcf_tbi}
		"""

rule gatk_filterMutect:
	input:
		vcf_tbi = OUTDIR + "/results/mutect2/{tumor}_vs_{patient}-N/unfiltered_{tumor}_vs_{patient}-N.vcf.idx",
		vcf = OUTDIR + "/results/mutect2/{tumor}_vs_{patient}-N/unfiltered_{tumor}_vs_{patient}-N.vcf",
		model = OUTDIR + "/results/mutect2/{tumor}_vs_{patient}-N/unfiltered_{tumor}_read_orientation_model.tar.gz",
		stats = OUTDIR + "/results/mutect2/{tumor}_vs_{patient}-N/unfiltered_{tumor}_vs_{patient}-N_merged.stats",
		interval = "/scratch/n/nicholsa/zyfniu/igenomes_ref/interval-files-folder/{num}-scattered.interval_list"
	output:
		filter_vcf = temp(OUTDIR + "/results/mutect2/{tumor}_vs_{patient}-N/filtered_{tumor}_vs_{patient}-N.{num}.vcf"),
		filter_vcf_index = temp(OUTDIR + "/results/mutect2/{tumor}_vs_{patient}-N/filtered_{tumor}_vs_{patient}-N.{num}.vcf.idx"),
		filter_vcf_stats = temp(OUTDIR + "/results/mutect2/{tumor}_vs_{patient}-N/filtered_{tumor}_vs_{patient}-N.{num}.vcf.filteringStats.tsv")
	threads: 2
	group: "variantCalling"
	shell:
		"""
		singularity exec -B $SCRATCH/igenomes_ref,{OUTDIR} \
			/gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/gatk-4.1.8.img \
		gatk --java-options "-Xmx8g" FilterMutectCalls \
		--ob-priors {input.model} \
		-stats {input.stats} \
		-R {REF_fasta} \
		-V {input.vcf} \
		-L {input.interval} \
		--output {output.filter_vcf}
		"""

def concat_vcf_filtered(wildcards):
	return	expand(OUTDIR + "/results/mutect2/" + wildcards.tumor + "_vs_" + wildcards.patient + "-N" + "/filtered_" + wildcards.tumor + "_vs_" + wildcards.patient + "-N.{num}.vcf", num = NUMBERS)

rule merge_mutect2_vcf_filtered:
	input:
		concat_vcf_filtered
	output:
		temp(OUTDIR + "/results/mutect2/{tumor}_vs_{patient}-N/filtered_{tumor}_vs_{patient}-N.vcf")
	threads: 2
	group: "variantCalling"
	shell:
		"singularity exec -B $SCRATCH/igenomes_ref,{OUTDIR} /gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/vcftools.img vcf-concat {input} > {output}"


rule index_filtered_vcf:
	input:
		vcf = OUTDIR + "/results/mutect2/{tumor}_vs_{patient}-N/filtered_{tumor}_vs_{patient}-N.vcf"
	output:
		vcf_tbi = OUTDIR + "/results/mutect2/{tumor}_vs_{patient}-N/filtered_{tumor}_vs_{patient}-N.vcf.idx"
	threads: 2
	group: "variantCalling"
	shell:
		"""
		singularity exec -B $SCRATCH/igenomes_ref,{OUTDIR} /gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/gatk-4.1.8.img \
		gatk --java-options "-Xmx8g" IndexFeatureFile \
		-I {input.vcf} \
		--output {output.vcf_tbi}
		"""

###
### VARIANT CALLING FOR ALL TUMOURS FROM THE SAME PATIENT USING COMBINED CALLING FROM MUTECT2
###

def get_bams_to_call(wildcards):
	samples = Multi[Multi.Patient == wildcards.patient].Sample
	return expand(OUTDIR +"/Recal/{tumor}.recal.bam", tumor = samples)

rule mutect2_all_tumours:
	input:
		interval = "/scratch/n/nicholsa/zyfniu/igenomes_ref/interval-files-folder/{num}-scattered.interval_list",
		normal = OUTDIR +"/Recal/{patient}-N.recal.bam",
		tumors = get_bams_to_call
	#	recurrence = "{sample}R/Recal/{sample}R.recal.bam"
	output:
		vcf = temp(OUTDIR + "/results/mutect2/all_{tumor}_vs_{patient}-N/unfiltered_{tumor}_vs_{patient}-N.{num}.vcf"),
		stats = temp(OUTDIR + "/results/mutect2/all_{tumor}_vs_{patient}-N/unfiltered_{tumor}_vs_{patient}-N.{num}.vcf.stats"),
		f12 = temp(OUTDIR + "/results/mutect2/all_{tumor}_vs_{patient}-N/unfiltered_{tumor}_vs_{patient}-N_f12.{num}.tar.gz"),
		index =  temp(OUTDIR + "/results/mutect2/all_{tumor}_vs_{patient}-N/unfiltered_{tumor}_vs_{patient}-N.{num}.vcf.idx")
	threads: 2
	group: "variantCalling"
	run:
		command = "singularity exec -B $SCRATCH/igenomes_ref,{OUTDIR} /gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/gatk-4.1.8.img \
		gatk --java-options \"-Xmx8g\" Mutect2 -R {REF_fasta} \
		-normal {wildcards.patient}-N \
		-L {input.interval} \
		--germline-resource {REF_gnomAD} \
		--panel-of-normals {REF_pon} --f1r2-tar-gz {output.f12} \
		-O {output.vcf}"
		for i in input[2:]: ###only the tumour files
			command = command + " -I " + i
		shell(command)

##SET NUMBERS FROM 0000 to 01999

NUMBERS = [str(a)+str(b)+str(c)+str(d) for a in range(0,1) for b in range(0,2) for c in range(0,10) for d in range(0,10)]

def concat_vcf_all_tumour(wildcards):
	return expand(OUTDIR + "/results/mutect2/all_" + wildcards.tumor + "_vs_" + wildcards.patient + "-N" + "/unfiltered_" + wildcards.tumor + "_vs_" + wildcards.patient + "-N.{num}.vcf", num = NUMBERS)

rule merge_mutect2_vcf_all_tumour:
	input:
		concat_vcf_all_tumour
	output:
		temp(OUTDIR + "/results/mutect2/all_{tumor}_vs_{patient}-N/unfiltered_{tumor}_vs_{patient}-N.vcf")
	threads: 1
	group: "variantCalling"
	shell:
		"""
		singularity exec -B $SCRATCH/igenomes_ref,{OUTDIR} /gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/vcftools.img \
		vcf-concat {input} > {output}
		"""

def concat_vcf_stats_all_tumour(wildcards):
	return expand(OUTDIR + "/results/mutect2/all_" + wildcards.tumor + "_vs_" + wildcards.patient + "-N" + "/unfiltered_" + wildcards.tumor + "_vs_" + wildcards.patient + "-N.{num}.vcf.stats", num = NUMBERS)

rule merge_stats_all_tumour:
	input:
		concat_vcf_stats_all_tumour
	output:
		temp(OUTDIR + "/results/mutect2/all_{tumor}_vs_{patient}-N/unfiltered_{tumor}_vs_{patient}-N_merged.stats")
	threads: 2
	group: "variantCalling"
	run:
		import os, glob
		cmd = "singularity exec -B $SCRATCH/igenomes_ref,{OUTDIR} /gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/gatk-4.1.8.img gatk MergeMutectStats "
		for input_file in input:
			cmd = cmd + " --stats " + input_file
		cmd = cmd + " -O {output}"
		shell(cmd)

def concat_vcf_f12_all_tumour(wildcards):
	return expand(OUTDIR + "/results/mutect2/all_" + wildcards.tumor + "_vs_" + wildcards.patient + "-N" + "/unfiltered_" + wildcards.tumor + "_vs_" + wildcards.patient + "-N_f12.{num}.tar.gz", num = NUMBERS)

rule gatk_LearnOrientationModel_all_tumour:
	input:
		concat_vcf_f12_all_tumour
	output:
		model = OUTDIR + "/results/mutect2/all_{tumor}_vs_{patient}-N/unfiltered_{tumor}_read_orientation_model.tar.gz"
	group: "variantCalling"
	threads: 2
	run:
		import os, glob
		cmd = "singularity exec -B $SCRATCH/igenomes_ref,{OUTDIR} /gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/gatk-4.1.8.img gatk LearnReadOrientationModel "
		for input_file in input:
			cmd = cmd + " -I " + input_file
		cmd = cmd + " -O {output}"
		shell(cmd)


rule index_unfiltered_vcf_all_tumour:
	input:
		vcf = OUTDIR + "/results/mutect2/all_{tumor}_vs_{patient}-N/unfiltered_{tumor}_vs_{patient}-N.vcf",
	output:
		vcf_tbi = OUTDIR + "/results/mutect2/all_{tumor}_vs_{patient}-N/unfiltered_{tumor}_vs_{patient}-N.vcf.idx"
	threads: 2
	group: "variantCalling"
	shell:
		"""
		singularity exec -B $SCRATCH/igenomes_ref,{OUTDIR} /gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/gatk-4.1.8.img \
		gatk --java-options "-Xmx8g" IndexFeatureFile \
		-I {input.vcf} \
		--output {output.vcf_tbi}
		"""

rule gatk_filterMutect_all_tumour:
	input:
		vcf_tbi = OUTDIR + "/results/mutect2/all_{tumor}_vs_{patient}-N/unfiltered_{tumor}_vs_{patient}-N.vcf.idx",
		vcf = OUTDIR + "/results/mutect2/all_{tumor}_vs_{patient}-N/unfiltered_{tumor}_vs_{patient}-N.vcf",
		model = OUTDIR + "/results/mutect2/all_{tumor}_vs_{patient}-N/unfiltered_{tumor}_read_orientation_model.tar.gz",
		stats = OUTDIR + "/results/mutect2/all_{tumor}_vs_{patient}-N/unfiltered_{tumor}_vs_{patient}-N_merged.stats",
		interval = "/scratch/n/nicholsa/zyfniu/igenomes_ref/interval-files-folder/{num}-scattered.interval_list"
	output:
		filter_vcf = temp(OUTDIR + "/results/mutect2/all_{tumor}_vs_{patient}-N/filtered_{tumor}_vs_{patient}-N.{num}.vcf"),
		filter_vcf_index = temp(OUTDIR + "/results/mutect2/all_{tumor}_vs_{patient}-N/filtered_{tumor}_vs_{patient}-N.{num}.vcf.idx"),
		filter_vcf_stats = temp(OUTDIR + "/results/mutect2/all_{tumor}_vs_{patient}-N/filtered_{tumor}_vs_{patient}-N.{num}.vcf.filteringStats.tsv")
	threads: 2
	group: "variantCalling"
	shell:
		"""
		singularity exec -B $SCRATCH/igenomes_ref,{OUTDIR} \
			/gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/gatk-4.1.8.img \
		gatk --java-options "-Xmx8g" FilterMutectCalls \
		--ob-priors {input.model} \
		-stats {input.stats} \
		-R {REF_fasta} \
		-V {input.vcf} \
		-L {input.interval} \
		--output {output.filter_vcf}
		"""

def concat_vcf_filtered_all_tumour(wildcards):
	return	expand(OUTDIR + "/results/mutect2/all_" + wildcards.tumor + "_vs_" + wildcards.patient + "-N" + "/filtered_" + wildcards.tumor + "_vs_" + wildcards.patient + "-N.{num}.vcf", num = NUMBERS)

rule merge_mutect2_vcf_filtered_all_tumour:
	input:
		concat_vcf_filtered_all_tumour
	output:
		temp(OUTDIR + "/results/mutect2/all_{tumor}_vs_{patient}-N/filtered_{tumor}_vs_{patient}-N.vcf")
	threads: 2
	group: "variantCalling"
	shell:
		"singularity exec -B $SCRATCH/igenomes_ref,{OUTDIR} /gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/vcftools.img vcf-concat {input} > {output}"


rule index_filtered_vcf_all_tumour:
	input:
		vcf = OUTDIR + "/results/mutect2/all_{tumor}_vs_{patient}-N/filtered_{tumor}_vs_{patient}-N.vcf"
	output:
		vcf_tbi = OUTDIR + "/results/mutect2/all_{tumor}_vs_{patient}-N/filtered_{tumor}_vs_{patient}-N.vcf.idx"
	threads: 2
	group: "variantCalling"
	shell:
		"""
		singularity exec -B $SCRATCH/igenomes_ref,{OUTDIR} /gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/gatk-4.1.8.img \
		gatk --java-options "-Xmx8g" IndexFeatureFile \
		-I {input.vcf} \
		--output {output.vcf_tbi}
		"""

rule annotate_mutect2_all_tumour:
	input:
		vcf = OUTDIR + "/results/mutect2/all_{tumor}_vs_{patient}-N/filtered_{tumor}_vs_{patient}-N.vcf"
	output:
		annotatedvcf = temp(OUTDIR + "/results/mutect2/all_{tumor}_vs_{patient}-N/{tumor}_vs_{patient}-N_snpEff.ann.vcf")
	threads: 2
	group: "variantCalling"
	shell:
		"""
		cd {OUTDIR}/results/mutect2/all_{wildcards.tumor}_vs_{wildcards.patient}-N
		singularity exec -B $SCRATCH/igenomes_ref,{OUTDIR} $SCRATCH/singularity_images/nfcore-sareksnpeff-2.6.GRCh38.img \
		snpEff -Xmx8g \
		GRCh38.86 \
		-csvStats {wildcards.tumor}_snpEff.csv \
		-nodownload \
		-canon \
		{input} \
		> {output}
#		mv snpEff_summary.html {wildcards.tumor}_snpEff.html
		"""

rule zip_snpeff_all_tumour:
	input:
		annotatedvcf = OUTDIR + "/results/mutect2/all_{tumor}_vs_{patient}-N/{tumor}_vs_{patient}-N_snpEff.ann.vcf"
	output: annotatedvcf = OUTDIR + "/results/mutect2/all_{tumor}_vs_{patient}-N/{tumor}_vs_{patient}-N_snpEff.ann.vcf.gz"
	threads: 2
	group: "variantCalling"
	shell:
		"""
		singularity exec -B $SCRATCH/igenomes_ref,{OUTDIR} \
			/gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/nfcore-sarek-2.6.img bgzip < {input} > {output}
		singularity exec -B $SCRATCH/igenomes_ref,{OUTDIR} \
			/gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/nfcore-sarek-2.6.img tabix {output}
		rm /scratch/n/nicholsa/zyfniu/AN_WGS/results/mutect2/all_{wildcards.tumor}_vs_{wildcards.patient}-N/*.idx -rf
		rm /scratch/n/nicholsa/zyfniu/AN_WGS/results/mutect2/all_{wildcards.tumor}_vs_{wildcards.patient}-N/*.filteringStats.tsv -rf
		"""

rule filter_mutect2_passOnly_all_tumour:
	input:
		annotatedvcf = OUTDIR + "/results/mutect2/all_{tumor}_vs_{patient}-N/{tumor}_vs_{patient}-N_snpEff.ann.vcf.gz"
	output: annotatedvcf = OUTDIR + "/results/mutect2/all_{tumor}_vs_{patient}-N/{tumor}_vs_{patient}-N_snpEff.ann.passOnly.vcf.gz"
	threads: 2
	group: "variantCalling"
	shell:
		"""
		singularity exec -B $SCRATCH/igenomes_ref,{OUTDIR} \
			/gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/nfcore-sarek-2.6.img bcftools view -f "PASS" {input} |\
		singularity exec -B $SCRATCH/igenomes_ref,{OUTDIR} \
		/gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/nfcore-sarek-2.6.img bgzip > {output}
		singularity exec -B $SCRATCH/igenomes_ref,{OUTDIR} \
		/gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/nfcore-sarek-2.6.img tabix -p vcf {output}
		"""


####
#### Manta
####
rule config_manta:
	input:
		normal = OUTDIR + "/Recal/{patient}-N.recal.bam",
		tumor = OUTDIR + "/Recal/{tumor}.recal.bam"
	output: OUTDIR + "/temp/Manta/{tumor}_vs_{patient}-N/runWorkflow.py"
	threads: 2
	group: "variantCalling"
	shell:
		"""
		singularity exec -B $SCRATCH/igenomes_ref,{OUTDIR} /gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/nfcore-sarek-2.6.img \
		configManta.py \
		--normalBam {input.normal} \
		--tumorBam  {input.tumor} \
		--reference {REF_fasta} \
		--runDir {OUTDIR}/temp/Manta/{wildcards.tumor}_vs_{wildcards.patient}-N
		"""

rule manta:
	input:
		normal = OUTDIR + "/Recal/{patient}-N.recal.bam",
		tumor = OUTDIR + "/Recal/{tumor}.recal.bam",
		script = OUTDIR + "/temp/Manta/{tumor}_vs_{patient}-N/runWorkflow.py"
	output:
		sv = OUTDIR + "/temp/Manta/{tumor}_vs_{patient}-N/results/variants/candidateSV.vcf.gz",
		smallindel = OUTDIR + "/temp/Manta/{tumor}_vs_{patient}-N/results/variants/candidateSmallIndels.vcf.gz",
		diploidSV = OUTDIR + "/temp/Manta/{tumor}_vs_{patient}-N/results/variants/diploidSV.vcf.gz",
		somaticSV = OUTDIR + "/temp/Manta/{tumor}_vs_{patient}-N/results/variants/somaticSV.vcf.gz",
		svIndex = OUTDIR + "/temp/Manta/{tumor}_vs_{patient}-N/results/variants/candidateSV.vcf.gz.tbi",
		smallindelIndex = OUTDIR + "/temp/Manta/{tumor}_vs_{patient}-N/results/variants/candidateSmallIndels.vcf.gz.tbi",
		diploidSVIndex = OUTDIR + "/temp/Manta/{tumor}_vs_{patient}-N/results/variants/diploidSV.vcf.gz.tbi",
		somaticSVIndex = OUTDIR + "/temp/Manta/{tumor}_vs_{patient}-N/results/variants/somaticSV.vcf.gz.tbi"
	group: "manta"
	threads: 80
	shell:
		"""
		singularity exec -B $SCRATCH/igenomes_ref,{OUTDIR} /gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/nfcore-sarek-2.6.img \
		python {input.script} -m local -j {threads} -g 160
		"""

rule mv_manta_files:
	input:
		file = OUTDIR + "/temp/Manta/{tumor}_vs_{patient}-N/results/variants/{structure}.vcf.gz",
		index = OUTDIR + "/temp/Manta/{tumor}_vs_{patient}-N/results/variants/{structure}.vcf.gz.tbi"
	output:
		file = OUTDIR + "/results/Manta/{tumor}_vs_{patient}-N/Manta_{tumor}_vs_{patient}-N.{structure}.vcf.gz",
		index = OUTDIR + "/results/Manta/{tumor}_vs_{patient}-N/Manta_{tumor}_vs_{patient}-N.{structure}.vcf.gz.tbi"
	group: "manta"
	threads: 1
	shell:
		"""
		cp {input.file} {output.file}
		cp {input.index} {output.index}
		"""
###
###	Annotatino using SNPEFF
###


rule annotate_mutect2:
	input:
		vcf = OUTDIR + "/results/mutect2/{tumor}_vs_{patient}-N/filtered_{tumor}_vs_{patient}-N.vcf"
	output:
		annotatedvcf = temp(OUTDIR + "/results/mutect2/{tumor}_vs_{patient}-N/{tumor}_vs_{patient}-N_snpEff.ann.vcf")
	threads: 2
	group: "variantCalling"
	shell:
		"""
		cd {OUTDIR}/results/mutect2/{wildcards.tumor}_vs_{wildcards.patient}-N
		singularity exec -B $SCRATCH/igenomes_ref,{OUTDIR} $SCRATCH/singularity_images/nfcore-sareksnpeff-2.6.GRCh38.img \
		snpEff -Xmx8g \
		GRCh38.86 \
		-csvStats {wildcards.tumor}_snpEff.csv \
		-nodownload \
		-canon \
		{input} \
		> {output}
#		mv snpEff_summary.html {wildcards.tumor}_snpEff.html
		"""

rule zip_snpeff:
	input:
		annotatedvcf = OUTDIR + "/results/mutect2/{tumor}_vs_{patient}-N/{tumor}_vs_{patient}-N_snpEff.ann.vcf"
	output: annotatedvcf = OUTDIR + "/results/mutect2/{tumor}_vs_{patient}-N/{tumor}_vs_{patient}-N_snpEff.ann.vcf.gz"
	threads: 2
	group: "variantCalling"
	shell:
		"""
		singularity exec -B $SCRATCH/igenomes_ref,{OUTDIR} \
			/gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/nfcore-sarek-2.6.img bgzip < {input} > {output}
		singularity exec -B $SCRATCH/igenomes_ref,{OUTDIR} \
			/gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/nfcore-sarek-2.6.img tabix {output}
		rm /scratch/n/nicholsa/zyfniu/AN_WGS/results/mutect2/{wildcards.tumor}_vs_{wildcards.patient}-N/*.idx -rf
		rm /scratch/n/nicholsa/zyfniu/AN_WGS/results/mutect2/{wildcards.tumor}_vs_{wildcards.patient}-N/*.filteringStats.tsv -rf
		"""

rule filter_mutect2_passOnly:
	input:
		annotatedvcf = OUTDIR + "/results/mutect2/{tumor}_vs_{patient}-N/{tumor}_vs_{patient}-N_snpEff.ann.vcf.gz"
	output: annotatedvcf = OUTDIR + "/results/mutect2/{tumor}_vs_{patient}-N/{tumor}_vs_{patient}-N_snpEff.ann.passOnly.vcf.gz"
	threads: 2
	group: "variantCalling"
	shell:
		"""
		singularity exec -B $SCRATCH/igenomes_ref,{OUTDIR} \
			/gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/nfcore-sarek-2.6.img bcftools view -f "PASS" {input} |\
		singularity exec -B $SCRATCH/igenomes_ref,{OUTDIR} \
		/gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/nfcore-sarek-2.6.img bgzip > {output}
		singularity exec -B $SCRATCH/igenomes_ref,{OUTDIR} \
		/gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/nfcore-sarek-2.6.img tabix -p vcf {output}
		"""

rule annotate_manta:
	input:
		vcf = OUTDIR + "/results/Manta/{tumor}_vs_{patient}-N/Manta_{tumor}_vs_{patient}-N.{structure}.vcf.gz"
	output:
		annotatedvcf = temp(OUTDIR + "/results/Manta/{tumor}_vs_{patient}-N/Manta_snpeff_{tumor}_vs_{patient}-N.{structure}.ann.vcf")
	threads: 2
	group: "manta"
	shell:
		"""
		cd {OUTDIR}/results/Manta/{wildcards.tumor}_vs_{wildcards.patient}-N
		singularity exec -B $SCRATCH/igenomes_ref,{OUTDIR} $SCRATCH/singularity_images/nfcore-sareksnpeff-2.6.GRCh38.img \
		snpEff -Xmx8g \
		GRCh38.86 \
		-csvStats {wildcards.tumor}_{wildcards.structure}_snpEff.csv \
		-nodownload \
		-canon \
		{input} \
		> {output}
#		mv snpEff_summary.html {wildcards.tumor}_{wildcards.structure}_snpEff.html
		"""

rule zip_manta:
	input:
		annotatedvcf = OUTDIR + "/results/Manta/{tumor}_vs_{patient}-N/Manta_snpeff_{tumor}_vs_{patient}-N.{structure}.ann.vcf"
	output: annotatedvcf = OUTDIR + "/results/Manta/{tumor}_vs_{patient}-N/Manta_snpeff_{tumor}_vs_{patient}-N.{structure}.ann.vcf.gz"
	threads: 2
	group: "manta"
	shell:
		"""
		singularity exec -B $SCRATCH/igenomes_ref,{OUTDIR} \
			/gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/nfcore-sarek-2.6.img bgzip < {input} > {output}
		singularity exec -B $SCRATCH/igenomes_ref,{OUTDIR} \
			/gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/nfcore-sarek-2.6.img tabix {output}
		"""

###
### COPY NUMBER USING ASCAT
###
rule alleleCount:
	input:
		bam = OUTDIR + "/Recal/{sample}.recal.bam",
		index = OUTDIR + "/Recal/{sample}.recal.bai",
		acloci = "/scratch/n/nicholsa/zyfniu/igenomes_ref/1000G_phase3_GRCh38_maf0.3.loci"
	output: OUTDIR + "/results/ASCAT/alleleCount/{sample}.alleleCount"
	threads: 2
	group: "variantCalling"
	shell:
		"""
		singularity exec -B $SCRATCH/igenomes_ref,{OUTDIR} /gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/nfcore-sarek-2.6.img \
		alleleCounter \
		-l {input.acloci} \
		-r {REF_fasta} \
		-b {input.bam} \
		-o {OUTDIR}/results/ASCAT/alleleCount/{wildcards.sample}.alleleCount
		"""

def getGender(wildcards):
		return expand("{gender}",gender = INPUT[INPUT.Patient == wildcards.patient].Sex.drop_duplicates())

###convert allele script from https://bitbucket.org/malinlarsson/somatic_wgs_pipeline/src/master/convertAlleleCounts.r
rule ConvertAlleleCounts:
	input:
		normal = OUTDIR + "/results/ASCAT/alleleCount/{patient}-N.alleleCount",
		tumor = OUTDIR + "/results/ASCAT/alleleCount/{tumor}.alleleCount",
		script = "/scratch/n/nicholsa/zyfniu/AN_WGS/AN_WGS_script/convertAlleleCounts.r"
	output:
		normalBaf = OUTDIR + "/results/ASCAT/{tumor}_vs_{patient}-N/{patient}-N.BAF",
		tumorBaf = OUTDIR + "/results/ASCAT/{tumor}_vs_{patient}-N/{tumor}.BAF",
		normalLogr = OUTDIR + "/results/ASCAT/{tumor}_vs_{patient}-N/{patient}-N.LogR",
		tumorLogr = OUTDIR + "/results/ASCAT/{tumor}_vs_{patient}-N/{tumor}.LogR"
	threads: 2
	params: gender = getGender
	group: "variantCalling"
	shell:
		"""
		cd {OUTDIR}/results/ASCAT/{wildcards.tumor}_vs_{wildcards.patient}-N
		singularity exec -B $SCRATCH/igenomes_ref,$SCRATCH/AN_WGS/AN_WGS_script,{OUTDIR} /gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/nfcore-sarek-2.6.img Rscript \
		/scratch/n/nicholsa/zyfniu/AN_WGS/AN_WGS_script/convertAlleleCounts.r \
		{wildcards.tumor} {input.tumor} \
		{wildcards.patient}-N {input.normal} \
		{params.gender}
		"""

rule ascat:
	input:
		normalBaf = OUTDIR + "/results/ASCAT/{tumor}_vs_{patient}-N/{patient}-N.BAF",
		tumorBaf = OUTDIR + "/results/ASCAT/{tumor}_vs_{patient}-N/{tumor}.BAF",
		normalLogr = OUTDIR + "/results/ASCAT/{tumor}_vs_{patient}-N/{patient}-N.LogR",
		tumorLogr = OUTDIR + "/results/ASCAT/{tumor}_vs_{patient}-N/{tumor}.LogR",
		acLociGC = "/scratch/n/nicholsa/zyfniu/igenomes_ref/1000G_phase3_GRCh38_maf0.3.loci.gc"
	output:
		results = OUTDIR + "/results/ASCAT/{tumor}_vs_{patient}-N/{tumor}_vs_{patient}-N.tumor.cnvs.txt"
	params:
		gender = getGender
	threads: 2
	group: "variantCalling"
	shell:
		"""
		for f in {OUTDIR}/results/ASCAT/{wildcards.tumor}_vs_{wildcards.patient}-N/*BAF \
		{OUTDIR}/results/ASCAT/{wildcards.tumor}_vs_{wildcards.patient}-N/*LogR; \
		do sed \'s/chr//g\' $f > {OUTDIR}/results/ASCAT/{wildcards.tumor}_vs_{wildcards.patient}-N/tmpFile; \
		mv {OUTDIR}/results/ASCAT/{wildcards.tumor}_vs_{wildcards.patient}-N/tmpFile $f;done
		cd {OUTDIR}/results/ASCAT/{wildcards.tumor}_vs_{wildcards.patient}-N
		singularity exec -B $SCRATCH/igenomes_ref,$SCRATCH/AN_WGS/AN_WGS_script,{OUTDIR} /gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/nfcore-sarek-2.6.img \
		Rscript /scratch/n/nicholsa/zyfniu/AN_WGS/AN_WGS_script/run_ascat.r \
        --tumorbaf {input.tumorBaf} \
        --tumorlogr {input.tumorLogr} \
        --normalbaf {input.normalBaf} \
        --normallogr {input.normalLogr} \
        --tumorname {wildcards.tumor} \
        --basedir {OUTDIR}/results/ASCAT/{wildcards.tumor}_vs_{wildcards.patient}-N \
        --gcfile {input.acLociGC} \
        --gender {params.gender}
		cd ../../../
		mv {OUTDIR}/results/ASCAT/{wildcards.tumor}_vs_{wildcards.patient}-N/{wildcards.tumor}.cnvs.txt {output.results}
		"""

###
###	END
###
