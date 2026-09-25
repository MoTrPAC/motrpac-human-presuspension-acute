# Figure 6 / Extended Data 8 / ST6 inputs

All from `MoTrPAC/precovid-analyses` PR #104 (branch `nmclark2/Figure7_v2.0`, commit
`e22945e`), `figures/landscape/figure_7/`, vendored verbatim. Every file takes the
`panel_source()` path; the two network tables also take an environment-variable override.

| File | What | Read by |
|---|---|---|
| `Gambardella-et-al-2020-table-s2.xlsx` | TFEB ChIP-seq peak annotation, Sci. Adv. 6, eabb0205, table S2 (sheet 2) | FIG6D, ED8A, ED8E |
| `TFEB_TARGET_GENES.v2024.1.Hs.tsv` | MSigDB C3 GTRD set TFEB_TARGET_GENES (Yevshin et al. 2019) | FIG6D, ED8A, ED8E |
| `GSM2354032_TFEB_e5_peaks.bed` | TFEB ChIP-seq peaks, hg19, Doronzo et al. 2019, EMBO J. 38, e98250 | FIG6D, ED8A, ED8E |
| `TFEB_targets_v2.xlsx` | the 705 SC-ION-predicted TFEB targets in the merged EE+RE network | FIG6D (ST6c), ED8A, FIG6E (check) |
| `SCION_supplemental_table_edges.csv` | every edge of the merged network after the 0.1 weight trim, 116,419 rows | ST6a, FIG6E; override `FIG6_SCION_EDGES_CSV` |
| `Merged_Network_combined_NMS_v2.csv` | Network Motif Score per node, 10,052 rows | ST6b; override `FIG6_SCION_NMS_CSV` |

The fourth ChIP source, Settembre et al. 2013 (Nat. Cell Biol. 15, 647), is two genes by
ChIP-qPCR (TFEB, PPARGC1A) and is coded in `figure_6/FIG6_ED8_helpers.R` rather than read from a file.

The SC-ION runs, the permutation trim and the Cytoscape merge that produced the two network
tables are not reproduced here. See [`docs/external_dependencies.md`](../../docs/external_dependencies.md).
