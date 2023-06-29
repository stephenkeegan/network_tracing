#!/usr/bin/env bash
### SLURM HEADER
##SBATCH --output=/projects/carter-lab/keegas/network_tracing/logs/batchtrace-%j.log  #!!
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=stephen.keegan@jax.org #!!
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
working="/projects/carter-lab/keegas/network_tracing"

cd $working
module load singularity

## exampl slurm batch job call:
#for i in results/biodomain/*/queryList*; do sbatch scripts/batch_trace.sbatch $i dir node edge; done

# Biodomain NWs
singularity exec $simg_path Rscript $working/scripts/Network_wKDA.R $1 $2 $3 $4 $5

# # Pseudotime NWs
# singularity exec $simg_path Rscript $working/scripts/pseudotime_nw_wKDA.R $1 $2 $3 $4 $5


exit
