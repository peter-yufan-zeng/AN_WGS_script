
#module purge
#module load singularity/3.8
#source ~/snakemake7.3.8/bin/activate

snakemake  --profile nichols_niagara \
 -s AN_WGS_script/run_wgs.snk \
 --configfile AN_WGS_script/niagara_config.yaml \
 --config input=$SCRATCH/20221219_gillison/20230112.csv outdir=$SCRATCH/20221219_gillison run_phylowgs=False \
 --singularity-args " -B /scratch/n/nicholsa/zyfniu,/scratch/n/nicholsa/zyfniu/gdc_reference,/scratch/n/nicholsa/zyfniu/20221219_gillison/reference " \
 --groups \
 sort_sam_to_bam=sort_sam_to_bam \
 merge_bam_mapped_and_index=merge_bam_mapped_and_index markdup=markdup \
 recalibrator=recalibrator gather_recal_tbl=recalibrator ApplyBQSR=recalibrator \
 Merge_Recal_Bam=Merge_Recal_Bam \
 samtools_stats=qc_samtools_stats bamqc=qc_bamqc \
 mutect2=mutect_single_variant_call mutect2_all_tumours=mutect_together_variant_call \
 --group-components \
 fastqc=12 sort_sam_to_bam=6 \
 merge_bam_mapped_and_index=6 markdup=3 \
 recalibrator=34 Merge_Recal_Bam=6 \
 qc_samtools_stats=20 qc_bamqc=20 \
 mutect_single_variant_call=20 mutect2_post_call=34 \
 mutect_together_variant_call=20 mutect2_all_post_call=34 \
 allelecount=10 ascat=10 \
 manta_call=8 manta_annotate=2 \
 --cores 80 \
 --resources mem_mb=160000 time=1440 threads=80 \
 --rerun-triggers input \
 --ri
