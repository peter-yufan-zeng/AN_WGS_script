def bwa_mem_fastq1(wildcards):
	return expand(INPUT[INPUT.Sample_Lane == wildcards.sample_lane].Fastq1)

def bwa_mem_fastq2(wildcards):
	return expand(INPUT[INPUT.Sample_Lane == wildcards.sample_lane].Fastq2)

rule hpviewer:
	input:
		fastq1 = bwa_mem_fastq1,
		fastq2 = bwa_mem_fastq2
	output:
		o1 = OUTDIR + "/results/HPViewer/{sample_lane}_HPV_profile.txt",
	threads: 80
	group: "hpv"
	shell:
		"""
		module load CCEnv
		module load StdEnv/2020
		module load scipy-stack/2020a bowtie2 samtools/1.10 bedtools
		source $SCRATCH/singularity_images/polyidus/bin/activate
		python $SCRATCH/AN_WGS/AN_WGS_script/scripts/HPViewer/HPViewer.py \
		-1 {input.fastq1} -2 {input.fastq2} -p {threads} -o {OUTDIR}/results/HPViewer/{wildcards.sample_lane}
		"""

rule polyidus:
    input:
		fastq1 = bwa_mem_fastq1,
		fastq2 = bwa_mem_fastq2
    output:
        o1 = OUTDIR + "results/{sample_lane}/polyidusOutput/results/HpvIntegrationInfo.tsv",
        i2 = OUTDIR + "results/{sample_lane}/polyidusOutput/results/exactHpvIntegrations.tsv"
    threads: 80
	group: "hpv"
	shell:
		"""
		module load CCEnv
		module load StdEnv/2020
		module load scipy-stack/2020a bowtie2 samtools/1.10 bedtools
		source $SCRATCH/singularity_images/polyidus/bin/activate
		python $SCRATCH/AN_WGS/AN_WGS_script/scripts/polyidus/polyidus.py $SCRATCH/igenomes_ref/bowtie2_index/Homo_sapiens_assembly38 \
		$SCRATCH/igenomes_ref/bowtie2_index/HPV16 \
		--fastq {input.r1} {input.r2} \
		--outdir {OUTDIR}/orphan/polyidus/{wildcards.sample}
		"""
