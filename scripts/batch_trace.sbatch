#!/usr/bin/env bash
### SLURM HEADER
#SBATCH --output=/projects/carter-lab/caryg/network_tracing/logs/batchtrace-%j.log  #!!
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=greg.cary@jax.org #!!
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
#simg_path="/projects/carter-lab/caryg/network_tracing/networktracing.sif"
simg_path="/projects/carter-lab/caryg/sif/rstudio_etc_4.2.1.sif"
working="/projects/carter-lab/caryg/network_tracing"

cd $working
module load singularity

#for i in {1..19}; do sbatch scripts/batch_trace.sbatch $i ; done
# singularity exec $simg_path Rscript $working/scripts/networktracing/code/05_Tracing_Pipeline.R $1

#singularity exec $simg_path Rscript $working/scripts/pathway_commons_tracing.R directed $1
#singularity exec $simg_path Rscript $working/scripts/pathway_commons_tracing.R undirected $1

singularity exec $simg_path Rscript $working/scripts/biodomain_network_tracing_and_kda.R $1 $2 $3 $4
#singularity exec $simg_path Rscript $working/scripts/biodomain_network_tracing_and_kda.R $1 dir
#singularity exec $simg_path Rscript $working/scripts/biodomain_network_tracing_and_kda.R $1 node
#singularity exec $simg_path Rscript $working/scripts/biodomain_network_tracing_and_kda.R $1 edge

exit
