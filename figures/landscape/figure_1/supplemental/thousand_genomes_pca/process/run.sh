#!/usr/bin/env bash

#SBATCH --job-name=process_vcf
#SBATCH --cpus-per-task=1
#SBATCH --partition=interactive
#SBATCH --account=default
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=1G

eval "$(micromamba shell hook --shell=bash)"
micromamba activate
micromamba activate snakemake

module load google-cloud-sdk
module load vcftools
module load plink2

BUCKET_PATH="gs://motrpac-portal-transfer-nygc/wgs/human/batch1_20220324/processed"

mkdir -p /projects/motrpac/PRECOVID/DATA/GENOTYPE/SAMPLES/

gsutil ls $BUCKET_PATH | grep Sample_[0-9]*/ > /projects/motrpac/PRECOVID/DATA/GENOTYPE/SAMPLES/sample_vcf_paths.txt

snakemake \
    --cluster "sbatch --cpus-per-task={threads} --job-name=process_vcf --mem-per-cpu={cluster.memory} --account={cluster.account} --time={cluster.time} --partition={cluster.partition}" \
    --jobs 8 \
    --cluster-config cluster_config.yaml \
    --rerun-trigger mtime
