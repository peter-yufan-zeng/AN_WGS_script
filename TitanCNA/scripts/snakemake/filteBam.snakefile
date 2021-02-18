configfile: "config/config.yaml"
configfile: "config/samples_filter.yaml"

rule all:
	input:
		expand("/scratch/n/nicholsa/zyfniu/singularity_images/TitanCNA/scripts/snakemake/results/{samples}.filtered.bam", samples=config["samples"])

rule filter_bam:
	input:
		lambda wildcards: config["samples"][wildcards.samples]
	output:
		"/scratch/n/nicholsa/zyfniu/singularity_images/TitanCNA/scripts/snakemake/results/{samples}.filtered.bam"
	params:
		chrs=config["chrs"]
	group:	"ichor"
	threads: 80
	shell:
		"""
		singularity exec -B \
		$SCRATCH/igenomes_ref,$SCRATCH/AN_WGS,\
$SCRATCH/singularity_images/TitanCNA/scripts/snakemake/results, \
 		/gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/nfcore-sarek-2.6.img samtools view \
		-L /scratch/n/nicholsa/zyfniu/igenomes_ref/hg38.chrom.bed -@ {threads} -o {output} {input}	
		singularity exec -B \
		$SCRATCH/igenomes_ref,$SCRATCH/AN_WGS,\
$SCRATCH/singularity_images/TitanCNA/scripts/snakemake/results, \
 		/gpfs/fs0/scratch/n/nicholsa/zyfniu/singularity_images/nfcore-sarek-2.6.img \
		samtools index {output}
		"""

