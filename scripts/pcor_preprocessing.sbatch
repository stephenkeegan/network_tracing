#!/usr/bin/env bash
### SLURM HEADER
#SBATCH --output=/projects/carter-lab/caryg/network_tracing/logs/pcor-%j.log  #!!
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=greg.cary@jax.org #!!
#SBATCH --qos=batch
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64GB
#SBATCH --export=ALL
### SLURM HEADER

localcores=${SLURM_CPUS_PER_TASK}

#simg_path="/home/caryg/rstudio_etc_4.0.3.sif"
simg_path="/projects/carter-lab/caryg/network_tracing/networktracing.sif"
working="/projects/carter-lab/caryg/network_tracing"

cd $working
module load singularity

singularity exec $simg_path Rscript $working/scripts/pcor_preprocessing.R

exit
