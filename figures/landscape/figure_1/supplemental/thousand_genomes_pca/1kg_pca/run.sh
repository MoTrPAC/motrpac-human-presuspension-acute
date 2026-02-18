#!/usr/bin/env bash

#SBATCH --job-name=1kg_pca
#SBATCH --cpus-per-task=1
#SBATCH --partition=interactive
#SBATCH --account=default
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=1G

eval "$(micromamba shell hook --shell=bash)"
micromamba activate
micromamba activate snakemake

module load tabix
module load bcftools

snakemake \
    --cluster "sbatch --cpus-per-task={threads} --job-name=1kg_pca --mem-per-cpu={cluster.memory} --account={cluster.account} --time={cluster.time} --partition={cluster.partition}" \
    --jobs 8 \
    --cluster-config cluster_config.yaml
