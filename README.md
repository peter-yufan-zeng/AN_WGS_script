# AN_WGS
The Anthony Nichols Lab WGS Pipeline

To run Compute Canada cluster, snakemake -s AN_WGS/Snakefile  --cores 1 -j 40 --cluster "sbatch -N 1 -t 20:00:00 --ntasks 80 --output=logs/%x_%j.log" --ri
