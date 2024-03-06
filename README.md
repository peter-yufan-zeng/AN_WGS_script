# AN_WGS: The Anthony Nichols Snakemake Pipelines for Next Generation Sequencing Data
The Anthony Nichols Lab pipeline, including scripts to run WGS, WES, RNASeq, scRNAseq, and more

To run Compute Canada Niagara cluster, snakemake -s AN_WGS/Snakefile  --cores 1 -j 40 --cluster "sbatch -N 1 -t 20:00:00 --ntasks 80 --output=logs/%x_%j.log" --ri

# To Set Up Snakemake on ComputeCanada Niagara Cluster

First, get access to the Niagara cluster and also [BBUFFER](https://docs.scinet.utoronto.ca/index.php/Burst_Buffer)

## Next, set up snakemake

```
module load python/3.8.5
cd ~
virtualenv create -n snakemake7.20.0
source ~/snakemake7.20.0/bin/activate
pip install pulp==2.7.0 snakemake==7.20.0 pandas
```
Copy the [profile config files](snakemake_configs/profiles/snakemake/) to your ~/.config/snakemake/ and make them executable with `chmod +x`

Set the following paths in your ~/.bashrc
```
export NXF_SINGULARITY_CACHEDIR=/path/to/singularity/cache
export AN_WGS_temp=$BBUFFER
export SNAKEMAKE_OUTPUT_CACHE=$SCRATCH/temp
export REF_DIR=/path/to/reference/directory 
```
##
Commands for individual pipelines could be found in [commands.txt](commands.txt)

# using mergerfs to create pooled disk space
'''
sudo mergerfs -o defaults,allow_other,use_ino,category.create=mfs,moveonenospc=true,minfreespace=1G /media/ionadmin/lab1:/media/ionadmin/lab2:/media/ionadmin/lab3:/media/ionadmin/lab4:/media/ionadmin/lab5:/media/ionadmin/lab6:/media/ionadmin/lab8:/media/ionadmin/lab9:/media/ionadmin/lab11:/media/ionadmin/lab13:/media/ionadmin/lab14:/media/ionadmin/lab15 /media/pool
'''

