#!/usr/bin/env bash
### SLURM HEADER
#SBATCH --output=/projects/carter-lab/caryg/network_tracing/logs/rstudio-%j.log  #!!
#SBATCH --mail-type=ALL
#SBATCH --mail-user=greg.cary@jax.org #!!
#SBATCH --qos=batch
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=50GB
#SBATCH --export=ALL
### SLURM HEADER

localcores=${SLURM_CPUS_PER_TASK}

#simg_path="/home/caryg/rstudio_etc_4.0.3.sif"
simg_path="/projects/carter-lab/caryg/network_tracing/networktracing.sif"
work_dir="/home/caryg/work"
d=$(pwd)

set -euo pipefail

#mkdir -p ${work_dir}
cd ${work_dir}
mkdir -p var var/run var/lib tmp

export PASSWORD=$(openssl rand -base64 15)
PORT=$(shuf -i8899-11999 -n1)

hostname_with_port=$(echo $(hostname -A):${PORT} | tr -d " ")

echo Login to ${hostname_with_port} with username: ${USER}  password: ${PASSWORD}

module load singularity

cd $d

singularity exec \
 -B $work_dir/var/run:/var/run/rstudio-server \
 -B $work_dir/var/lib:/var/lib/rstudio-server \
 -B $work_dir/tmp:/tmp \
 -B /projects/carter-lab/caryg/ \
 ${simg_path} rserver \
 --www-port ${PORT} \
 --auth-none=0 --auth-pam-helper-path=pam-helper \
 --auth-timeout-minutes=0 \
 --auth-stay-signed-in-days=30
