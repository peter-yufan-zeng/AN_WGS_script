snakemake  --profile nichols_niagara \
 -s AN_WGS_script/snakemake_scripts/scrnaseq_velocyto.snk \
 --configfile AN_WGS_script/niagara_config_pez.yaml \
 --config input=$SCRATCH/20230227_isgs_scrnaseq/input.csv outdir=$SCRATCH/20230227_isgs_scrnaseq \
 --groups velocyto=cellranger \
 --group-components \
 cellranger=12 \
 --rerun-triggers input \
 --ri
