#!/bin/bash
# Run the pre-COVID landscape cross-tissue PLIER RNA analysis

# Exit immediately if a command exits with a non-zero status
set -e

# Run the R script
Rscript precovid_landscape_Figure5_crosstissueplierrna.R

