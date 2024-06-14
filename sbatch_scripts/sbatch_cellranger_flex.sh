#!/bin/bash
#
#
# =============================================================================
# Cellranger_Flex Script
# =============================================================================
#
## Run with the command sbatch --time=$(awk 'END { printf "%02d:00:00", int((NR - 2) * 0.9 + 0.99) }' "Input.csv") Cellranger_Flex.sh
# Settings and Parameters
###this script uses cellranger 8.0.1 and refdata-gex-GRCh38-2020-A
Organism=$(head -n 1 "Input.csv" | tr -d '[:space:]')
FASTQs=$(sed -n '2p' Input.csv | tr -d '[:space:]')
if [ "$Organism" == "human" ];	then
	echo "[gene-expression]
reference,$REF_DIR/cellranger/refdata-gex-GRCh38-2020-A
probe-set,$REF_DIR/cellranger/Chromium_Human_Transcriptome_Probe_Set_v1.0_GRCh38-2020-A.csv
create-bam,true

[libraries]
fastq_id,fastqs,feature_types
Pool_44,$FASTQs,Gene Expression

[samples]
sample_id,probe_barcode_ids,description" > config.csv
	tail -n +3 "Input.csv" >> config.csv
elif [ "$Organism" == "mouse" ];   then
        echo "[gene-expression]
reference,$REF_DIR/cellranger/refdata-gex-mm10-2020-A
probe-set,$REF_DIR/cellranger/Chromium_Mouse_Transcriptome_Probe_Set_v1.0_mm10-2020-A.csv
create-bam,true

[libraries]
fastq_id,fastqs,feature_types
Pool_44,$FASTQs,Gene Expression

[samples]
sample_id,probe_barcode_ids,description" > config.csv
	tail -n +3 "Input.csv" >> config.csv
else
	echo "invalid input"
	exit 1
fi

# Directives
#SBATCH --job-name Cellranger_Flex
#SBATCH --export=ALL
#SBATCH --nodes=1 --ntasks-per-node=40
#SBATCH --signal=2
#SBATCH --no-requeue

$REF_DIR/cellranger/cellranger-8.0.0/cellranger multi --id=cellranger_outputs --csv=config.csv --jobmode=local --localcores=40 --localmem=170
