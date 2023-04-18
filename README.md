# AN_WGS: The Anthony Nichols WGS Pipeline
The Anthony Nichols Lab WGS Pipeline

To run Compute Canada Niagara cluster, snakemake -s AN_WGS/Snakefile  --cores 1 -j 40 --cluster "sbatch -N 1 -t 20:00:00 --ntasks 80 --output=logs/%x_%j.log" --ri
# using mergerfs to create pooled disk space
sudo mergerfs -o defaults,allow_other,use_ino,category.create=mfs,moveonenospc=true,minfreespace=1G /media/ionadmin/lab1:/media/ionadmin/lab2:/media/ionadmin/lab3:/media/ionadmin/lab4:/media/ionadmin/lab5:/media/ionadmin/lab6:/media/ionadmin/lab8:/media/ionadmin/lab9:/media/ionadmin/lab11:/media/ionadmin/lab12:/media/ionadmin/lab13:/media/ionadmin/lab14:/media/ionadmin/lab15 /media/pool
