#!/bin/bash

# Ensure we receive the parameters
PERMUTE=$1
NUM_CORES=$2

set -e
module load R/4.2.2

# Pass the parameters to the R Markdown rendering process
Rscript -e "rmarkdown::render('scion_figures.Rmd', params=list(permute=${PERMUTE}, num_cores=${NUM_CORES}))"

