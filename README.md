# MotrpacPreSuspensionAcute

This repository tracks the exact code, parameters, and documentation/links to external data used to generate each of the manuscripts for the Molecular Transducers of Physical Activity Consortium Pre-Suspension Human Phase.

For each manuscript, the code used to generate the figures is located in the subfolder `figures/`, labeled under the relevant manuscript name. The `landscape` folder indicates the integrative analysis that incorporates data from all tissues. Each topic sub-analysis group is responsible for the relevant code for their figures. If you have any questions or issues that show up, please [submit a new issue](https://github.com/MoTrPAC/MotrpacPreSuspensionAcute/issues){target="_blank"} 
and be very clear about which specific code or figure panel you are describing. Please include as many details as possible. 

This repository also contains some of the quality control steps that visualize decisions such as how outliers were flagged, or how the decisions were made that reflect the data generation process that is implemented in `MotrpacHumanPreSuspensionAnalysis`.
Please review the`QC` folder to take a look at these visualizations and data processing steps, which will explain in a bit more detail. 


## Repository Structure

```
MotrpacPreSuspensionAcute/

├── figures/                           # All manuscript figure code and outputs
│   ├── landscape/                     # Cross-tissue integrative analysis paper
│   ├── blood/                         # Blood tissue paper
│   ├── adipose/                       # Adipose tissue paper
│   ├── muscle/                        # Muscle tissue paper
│   └── splicing/                      # Alternative splicing paper
└── QC/                                # Quality control vignettes (HTML reports)
```

## Notes for the public:

If you are viewing the original bioRxiv submission, please realize that the the code used to generate the figures may change over time as the consortium analysis groups address revisions or concerns or do extra analysis. Please refer to the `RELEASES` list on the side of the github repository for more information on if you would like stored versions of the repository corresponding to specific versions of the bioRxiv submission. As of Feb 25th, 2026, the code used for each of the `landscape`, `splicing`, and `adipose` figures is available, but analysts are still updating the code used for the blood and muscle manuscripts. The first formal github release will be locked in once all the code for the panels for those manuscripts is available. 

To protect participant privacy and comply with data-use governance policies, individual-level (subject-level) molecular or phenotypic data, found in `MotrpacHumanPreSuspensionData`, are available only through formal data access requests to the MoTrPAC consortium. This means that some of the figures generated will not be directly replicate-able until that access is granted.

The figures folders rely HEAVILY on the the packages for `MotrpacHumanPreSuspensionAnalysis` and `MotrpacHumanPreSuspensionData` for access. Some references to data from external datasets are located within the relevant figure folder, but most of the data loaded will be done through the `MotrpacHumanPreSuspensionData` package (which is available via request to the Consortium). 



