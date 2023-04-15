#!/usr/bin/env bash
### SLURM HEADER
##SBATCH --output=/projects/carter-lab/caryg/network_tracing/logs/batchtrace-%j.log  #!!
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=greg.cary@jax.org #!!
##SBATCH --qos=batch
##SBATCH --qos=long
##SBATCH --time=6-00:00:00
##SBATCH --nodes=1
##SBATCH --ntasks=1
##SBATCH --cpus-per-task=12
##SBATCH --mem=16GB
##SBATCH --export=ALL
### SLURM HEADER

localcores=${SLURM_CPUS_PER_TASK}

#simg_path="/home/caryg/rstudio_etc_4.0.3.sif"
#simg_path="/projects/carter-lab/caryg/network_tracing/networktracing.sif"
simg_path="/projects/carter-lab/caryg/sif/rstudio_etc_4.2.1.sif"
working="/projects/carter-lab/caryg/network_tracing"

cd $working
module load singularity

## exampl slurm batch job call:
#for i in {1..19}; do sbatch scripts/batch_trace.sbatch $i dir node edge; done

#singularity exec $simg_path Rscript $working/scripts/biodom_nw_path_trace.R $1 $2 $3 $4
singularity exec $simg_path Rscript $working/scripts/biodom_nw_wKDA.R $1 $2 $3 $4 $5


exit
