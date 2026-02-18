# VCF File Processing

## Processing Steps

1. Enter directory using `cd process/` and run `sbatch run.sh` to process WGS files.
2. Enter directory using `cd 1kg_pca/` and run `sbatch run.sh` to create PCA of genotypes with the 1000 Genomes project.
3. Run through `vcf_quality_metrics.Rmd` to generate QC plots for the genotypes.
4. Run through `1kg_pca.Rmd` to generate PCA plots for the genotypes.

