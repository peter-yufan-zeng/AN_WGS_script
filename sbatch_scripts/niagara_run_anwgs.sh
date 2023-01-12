
module purge
module load singularity/3.8
source ~/snakemake7.3.8/bin/activate

snakemake  --profile nichols_niagara \
	-s AN_WGS_script/run_wgs.snk \
	--config input=$SCRATCH/20221219_gillison/20230110_fastq.csv outdir=$SCRATCH/20221219_gillison \
	--singularity-args " -B $SCRATCH/gdc_reference " \
  --group-components \ # this groups jobs together on the same cluster to save
  ### wgs_align
  fastqc=12 bwa_mem=3 sort_sam_to_bam=6 \
	### wgs_markdup_recal
	merge_bam_mapped_and_index=6 markdup=3 \
	recalibrator=34 Merge_Recal_Bam=6 \
	###bam_qc
	qc_samtools_stats=20 qc_bamqc=20 \
	### wgs_mutect2_individual
	mutect_single_variant_call=20 mutect2_post_call=34 \
	### wgs_mutect2_together
	mutect_together_variant_call=20 mutect2_all_post_call=34 \
	### wgs_ascat
  allelecount=10 ascat=10 \
	### wgs_manta
	manta_call=8 manta_annotate=2 \
	--max-threads 80 \
	--ri
