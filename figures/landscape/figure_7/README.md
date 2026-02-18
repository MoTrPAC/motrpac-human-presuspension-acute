# Overview

This folder contains all files required to reproduce the SCION network analysis and figure panels included in the manuscript. The SCION workflow integrates transcriptomic and proteomic data to infer regulatory networks and identify condition-specific signaling relationships. The files here allow full regeneration of both the raw analytical outputs (via the SCION algorithm) and the final publication-quality figure panels.

## Reproducibility Workflow

### Primary Analysis and Data Preparation — scion_figures.Rmd

This R Markdown file serves as the main analysis script. It:

Loads and preprocesses the raw SCION input data (e.g., gene expression matrices, metadata, prior network structures).

Executes the SCION algorithm to infer regulatory and signaling connections across tissues and experimental conditions.

Saves key intermediate objects (e.g., inferred edges, node-level attributes, confidence scores, and pathway annotations) for downstream visualization.

Output files from this step are used as input for Cytoscape and the figure generation scripts below.

Running this document ensures that all computations are reproducible from the original data sources.

### Network Visualization and Quality Review — Cytoscape

After running scion_figures.Rmd, import the resulting edge and node tables into Cytoscape.

Use Cytoscape (see: `Figure7_SCION.ai`) to:

see how network identification was formally preformed

### Final Figure Assembly — Figure_7_script.R and Figure_S7_script.R

These scripts recreate the final figure panels used in the manuscript:

Figure_7_script.R generates the main text SCION figure, integrating Cytoscape layouts with supplementary statistical overlays and annotations.

Figure_S7_script.R generates the supplementary figure panels, including expanded network contexts, node attribute plots, and comparative analyses.

Both scripts pull in the processed SCION output (from Step 1) and any exported layout coordinates or attributes from Cytoscape.

These figures are then contextualized against relevant literature sources to highlight biological significance and known signaling interactions.
