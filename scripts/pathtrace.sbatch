#!/usr/bin/env bash
### SLURM HEADER
##SBATCH --output=/projects/carter-lab/caryg/network_tracing/logs/batchtrace-%j.log  #!!
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=greg.cary@jax.org #!!
#SBATCH --qos=batch
##SBATCH --time=96:00:00
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
#for i in results/biodomain/*/queryList*; do sbatch scripts/batch_trace.sbatch $i dir node edge; done

# Biodomain NWs
singularity exec $simg_path Rscript $working/scripts/Network_path_trace.R $1 $2 $3 $4 $5

# Biodomain NWs
singularity exec $simg_path Rscript $working/scripts/biodom_filter_traced_nw.R $1 $2 $3 $4 $5

# # Pseudotime NWs
# singularity exec $simg_path Rscript $working/scripts/pseudotime_nw_path_trace.R $1 $2 $3 $4 $5

exit
