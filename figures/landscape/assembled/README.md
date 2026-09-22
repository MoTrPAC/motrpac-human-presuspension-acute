# Assembled figures and supplementary tables

The figures and supplementary tables as submitted. Everything here is an **output**, kept
for reference: nothing in this folder is read by any script, and nothing regenerates it.

## How these relate to the scripts

Each panel in these PDFs was drawn by the figure script named below, which writes one PDF per
panel to `outputs/panels/<figure_dir>/`. The assembled file is those panels placed together on
a page.

**The assembled panels are not always byte-identical to what the scripts emit.** Placing the
panels is done in Illustrator, and small text was adjusted there: font sizes, label positions,
a truncated axis label given room, panel letters added. **No values, no statistics and no
graphical elements were changed**: a dot in an assembled figure sits where the script put it.
So a panel rebuilt from the script may differ cosmetically from its counterpart here, and that
difference is typesetting rather than analysis.

Five panels are not drawn by any script and exist only in the assembled files:

| Panel | Made with |
|---|---|
| FIG1A | BioRender |
| FIG1B | BioRender |
| FIG6A | Illustrator |
| FIG6B | Cytoscape |
| FIG6G | BioRender |

## Figures

| Assembled file | Built by | Panels |
|---|---|---|
| `Figure1_Overview.pdf` | `figure_1/FIG1.R` | FIG1C |
| `ExtDataFig1_Overview_supp.pdf` | `figure_1/ED1.R` | ED1A, ED1B, ED1C, ED1D |
| `Figure2_DA.pdf` | `figure_2/FIG2.R` | FIG2A, FIG2B, FIG2Bii, FIG2C, FIG2D, FIG2E, FIG2F |
| `ExtDataFig2_DA_supp.pdf` | `figure_2/ED2.R` | ED2A, ED2B, ED2C, ED2D, ED2E, ED2F |
| `ExtDataFig3_SexEffects.pdf` | `extended_data_3/ED3.R` | ED3A, ED3B, ED3C, ED3D, ED3E |
| `Figure3_RE_vs_EE.pdf` | `figure_3/FIG3.R` | FIG3A, FIG3B, FIG3C, FIG3D, FIG3EFG |
| `ExtDataFig4_RE_vs_EE_supp.pdf` | `figure_3/ED4.R` | ED4A, ED4B, ED4C, ED4D, ED4E, ED4F, ED4G |
| `Figure4_cmeans_PLIER.pdf` | `figure_4/FIG4.R` | FIG4A, FIG4B, FIG4C, FIG4D, FIG4E, FIG4F, FIG4G, FIG4H, FIG4I |
| `ExtDataFig5_clustering_supp.pdf` | `figure_4/ED5.R` | ED5A, ED5B, ED5C, ED5D, ED5E |
| `ExtDataFig6_PLIER_supp.pdf` | `figure_4/ED6.R` | ED6A, ED6Aii, ED6B, ED6C, ED6D, ED6E, ED6F, ED6G, ED6H, ED6I |
| `Figure5_HOMER.pdf` | `figure_5/FIG5.R` | FIG5A, FIG5B, FIG5C, FIG5D, FIG5E, FIG5F, FIG5G, FIG5H |
| `ExtDataFig7_HOMER_supp.pdf` | `figure_5/ED7.R` | ED7A, ED7B |
| `Figure6_SCION.pdf` | `figure_6/FIG6.R` | FIG6C, FIG6D, FIG6E, FIG6F |
| `ExtDataFig8_SCION_supp.pdf` | `figure_6/ED8.R` | ED8A, ED8B, ED8C, ED8D, ED8E, ED8F |
| `Figure7_Secretome.pdf` | `figure_7/FIG7.R` | FIG7A, FIG7B, FIG7C, FIG7D |
| `ExtDataFig9_Secretome_supp.pdf` | `figure_7/ED9.R` | ED9A, ED9B, ED9C, ED9D |
| `SuppFigS1_ATACseqQC.pdf` | `supplementary_figure_1/SF1.R` | SF1A, SF1B, SF1C, SF1D, SF1E, SF1F, SF1G |

## Supplementary tables

Each workbook is one supplementary table number, with a sheet per sub-table. The scripts write
one tab-delimited file per sub-table to `outputs/tables/`; assembling those into workbooks is a
submission step and is not reproduced here.

| Workbook | Sub-tables | Note |
|---|---|---|
| `Landscape_Supplementary_Table_1.xlsx` | ST1a, ST1b, ST1c, ST1d, ST1e | all built here |
| `Landscape_Supplementary_Table_2.xlsx` | ST2a, ST2b, ST2c, ST2d, ST2e, ST2f, ST2g, ST2h | all built here |
| `Landscape_Supplementary_Table_3.xlsx` | ST3a, ST3b, ST3c, ST3d, ST3e, ST3f, ST3g, ST3h | ST3e not built here |
| `Landscape_Supplementary_Table_4.xlsx` | ST4a, ST4b, ST4c, ST4d, ST4e | ST4b-ST4e not built here |
| `Landscape_Supplementary_Table_5.xlsx` | ST5a | all built here |
| `Landscape_Supplementary_Table_6.xlsx` | ST6a, ST6b, ST6c, ST6d, ST6e, ST6f, ST6g | all built here |
| `Landscape_Supplementary_Table_7.xlsx` | ST7a | all built here |

Extended Data Table 1 is in the manuscript body rather than the supplement, and has no workbook
here. `figure_3/EDT1.R` writes it.

`config/panel_map.json` and `config/table_map.json` are authoritative for which script owns
which panel and which sub-table.
