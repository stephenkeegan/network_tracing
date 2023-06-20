#!/usr/bin/env bash
### SLURM HEADER
#SBATCH --output=/projects/carter-lab/caryg/network_tracing/logs/master_batchtrace-%j.log  #!!
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=greg.cary@jax.org #!!
#SBATCH --qos=batch
#SBATCH --export=ALL
### SLURM HEADER

working="/projects/carter-lab/caryg/network_tracing"

cd $working

## exampl slurm batch job call:
#for i in results/biodomain/*/queryList*; do sbatch scripts/batch_trace.sbatch $i dir node edge; done

## wait for trace to complete
trace_jobid=`sbatch --job-name=pathtrace -o ./logs/trace-%j.log --time=72:00:00 --cpus-per-task=12 --mem=16GB ./scripts/pathtrace.sbatch $1 $2 $3 $4 | grep -o "[0-9]*"`
echo trace_jobid

kda_jobid_ef1=`sbatch --job-name=kda -o ./logs/kda-%j.log --dependency=afterok:$trace_jobid --qos=long --time=14-00:00:00 --cpus-per-task=1 --mem=4GB ./scripts/kda.sbatch $1 $2 $3 $4 ef1`
echo kda_jobid_ef1

kda_jobid_ef2=`sbatch --job-name=kda -o ./logs/kda-%j.log --dependency=afterok:$trace_jobid --qos=long --time=14-00:00:00 --cpus-per-task=1 --mem=4GB ./scripts/kda.sbatch $1 $2 $3 $4 ef0`
echo kda_jobid_ef2

## don't wait, just run kda:
# kda_jobid_ef1=`sbatch --job-name=kda -o ./logs/kda-%j.log --qos=long --time=14-00:00:00 --cpus-per-task=1 --mem=4GB ./scripts/kda.sbatch $1 $2 $3 $4 ef1`
# echo kda_jobid_ef1
# 
# kda_jobid_ef0=`sbatch --job-name=kda -o ./logs/kda-%j.log --qos=long --time=14-00:00:00 --cpus-per-task=1 --mem=4GB ./scripts/kda.sbatch $1 $2 $3 $4 ef0`
# echo kda_jobid_ef0

exit
