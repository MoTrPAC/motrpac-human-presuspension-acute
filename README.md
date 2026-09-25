# motrpac-human-presuspension-acute

This repository tracks the exact code, parameters, and documentation and links to external data used to generate each manuscript from the Molecular Transducers of Physical Activity Consortium (MoTrPAC) Pre-Suspension human phase.

For each manuscript, the code used to generate the figures is located in the subfolder `figures/`, labeled under the relevant manuscript name. The `landscape` folder contains the integrative analysis that incorporates data from all tissues. Each topic sub-analysis group is responsible for the relevant code for their figures. If you have questions or find a problem, please [submit a new issue](https://github.com/MoTrPAC/motrpac-human-presuspension-acute/issues) and state which script or figure panel you are describing, with as much detail as possible.

This repository also contains quality control reports that document decisions such as how outliers were flagged and how the data processing implemented in `MotrpacHumanPreSuspensionAnalysis` was chosen. See the `QC` folder for these reports.

## How this repo fits with the others

The Pre-Suspension human work is split across four repositories. Individual-level molecular
and phenotypic data cannot be distributed publicly, so they live in a separate, access-gated
package, and everything that *can* be released publicly (aggregate results, all analysis
code) is kept clear of them.

```
                  MoTrPAC BIC — consortium GCS buckets
                        (raw assay data, gated)
                                   │
                   motrpac-human-presuspension-repro
           normalizes omics data, applies statistical models,
           builds every data object, versions it, uploads it,
                   and carries it into both packages
                                   │
                ┌──────────────────┴──────────────────┐
                ▼                                     ▼
  MotrpacHumanPreSuspensionData       MotrpacHumanPreSuspensionAnalysis
subject-level data — access-gated         aggregate results — public
                └──────────────────┬──────────────────┘
                                   ▼
                   motrpac-human-presuspension-acute
                 manuscript figure code + QC vignettes
                                   ▼
                              manuscripts
```

| Repository | What it holds | Access |
|---|---|---|
| [`motrpac-human-presuspension-repro`](https://github.com/MoTrPAC/motrpac-human-presuspension-repro) | the end-to-end rebuild pipeline and its pinned software environment | code; a full run needs consortium bucket access |
| [`MotrpacHumanPreSuspensionData`](https://github.com/MoTrPAC/MotrpacHumanPreSuspensionData) | subject-level molecular and phenotypic data objects | formal data-access request to the consortium |
| [`MotrpacHumanPreSuspensionAnalysis`](https://github.com/MoTrPAC/MotrpacHumanPreSuspensionAnalysis) | differential analysis, group summary statistics, enrichment, clustering, feature-to-gene map, plotting functions | public |
| [`motrpac-human-presuspension-acute`](https://github.com/MoTrPAC/motrpac-human-presuspension-acute) | per-manuscript figure code and QC vignettes | code public; some panels need `MotrpacHumanPreSuspensionData` access |

## Repository Structure

```
motrpac-human-presuspension-acute/
├── figures/                           # All manuscript figure code and outputs
│   ├── landscape/                     # Cross-tissue integrative analysis paper
│   ├── blood/                         # Blood tissue paper
│   ├── adipose/                       # Adipose tissue paper
│   ├── muscle/                        # Muscle tissue paper
│   └── splicing/                      # Alternative splicing paper
└── QC/                                # Quality control reports (HTML)
```

## Notes for the public

The code used to generate the figures may change over time as the consortium analysis groups address revisions, respond to reviewer concerns, or add analyses. The version of this repository that matches the original bioRxiv submission is tagged as the [`v1.0-alpha` release](https://github.com/MoTrPAC/motrpac-human-presuspension-acute/releases/tag/v1.0-alpha); see the Releases list in the GitHub sidebar for stored versions. Code for all five manuscripts (`landscape`, `blood`, `adipose`, `muscle` and `splicing`) is in `figures/`. The tissue papers (`adipose`, `blood` and `muscle`) are still a work in progress, and their code may change.

Since that release, the `landscape` code on `main` has been reorganized around the manuscript's own figure numbering: one script per figure, run together by `figures/landscape/run_all.R`, and aligned with `MotrpacHumanPreSuspensionAnalysis` 2.0.8 and `MotrpacHumanPreSuspensionData` 2.0.3. [`figures/landscape/README.md`](figures/landscape/README.md) gives the script for every figure panel and supplementary table.

To protect participant privacy and comply with data-use governance policies, individual-level (subject-level) molecular and phenotypic data, found in `MotrpacHumanPreSuspensionData`, are available only through formal data access requests to the MoTrPAC consortium. Some figures therefore cannot be reproduced until that access is granted.

The figure folders rely heavily on the `MotrpacHumanPreSuspensionAnalysis` and `MotrpacHumanPreSuspensionData` packages. Some data from external datasets are stored in the relevant figure folder, but most data are loaded through `MotrpacHumanPreSuspensionData`.

## Running the figure scripts

Instructions for running each manuscript's figure code:

- `landscape`: [`figures/landscape/README.md`](figures/landscape/README.md)
- `adipose`: [`figures/adipose/README.md`](figures/adipose/README.md) (work in progress)
- `muscle`: [`figures/muscle/README.md`](figures/muscle/README.md) (work in progress)
- `blood`: [`figures/blood/README.md`](figures/blood/README.md) (work in progress)
- `splicing`: [`figures/splicing/`](figures/splicing/)
