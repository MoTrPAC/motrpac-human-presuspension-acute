# PreSuspension-Analyses


Repository for the collaborative multi-omic analysis of Pre-Covid data as part of MoTrPAC. 

This repository contains the exact code and parameters used to regenerate the figures included in the manuscript, located in the `/figures/` folder.

This repository also contains some of the quality control steps that visualize decisions such as how outliers were flagged, or how the decisions were made that reflect the data generation process that is implemented in `MotrpacHumanPreSuspensionAnalysis`. 

## Reviewer/Public Note:

To protect participant privacy and comply with data-use governance policies, individual-level (subject-level) molecular or phenotypic data are available only through formal data access requests to the MoTrPAC consortium. This means that some of the figures generated will not be directly replicate-able until that access is granted. 

For any public member evaluating the code for any of the MoTrPAC PreSuspension phase:
please refer to the README files found in `figures` folder. Each subfolder in the figures folder will correspond to the analysis for the relevant manuscript. The 'landscape' folder corresponds to the manuscript for the cross-tissue integrative analysis, whereas each tissue-focused paper is named after the relevant tissue. 

The figures folders rely HEAVILY on the the packages for `MotrpacHumanPreSuspensionAnalysis` and `MotrpacHumanPreSuspensionData` for access. Some references to data from external datasets are located within the relevant figure folder, but most of the data loaded will be done through the `MotrpacHumanPreSuspensionData` package (which is available via request). 

The README in each `/figures/*` folder should describe if each individual panel requires individual-level access.


