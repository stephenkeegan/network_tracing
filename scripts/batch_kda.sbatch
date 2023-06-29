#!/usr/bin/env bash
### SLURM HEADER
#SBATCH --output=/projects/carter-lab/keegas/network_tracing/logs/batchKDA-%j.log  #!!
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=stephen.keegan@jax.org #!!
#SBATCH --qos=batch
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64GB
#SBATCH --export=ALL
### SLURM HEADER

localcores=${SLURM_CPUS_PER_TASK}

#simg_path="/home/caryg/rstudio_etc_4.0.3.sif"
simg_path="/projects/carter-lab/caryg/network_tracing/networktracing.sif"
working="/projects/carter-lab/keegas/network_tracing"

cd $working
module load singularity

#for i in {1..16}; do sbatch scripts/batch_kda.sbatch $i ; done
singularity exec $simg_path Rscript $working/scripts/wKDA.R $1

exit
