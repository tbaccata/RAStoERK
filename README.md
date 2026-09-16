[![GPLv3 license](https://img.shields.io/badge/License-GPLv3-blue.svg)](http://perso.crans.org/besson/LICENSE.html)
[![bioRxiv](https://img.shields.io/badge/bioRxiv-10.64898/2026.04.07.716061v1-b31b1b.svg)](https://biorxiv.org/content/10.64898/2026.04.07.716061v1)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22798879.svg)](https://doi.org/10.5281/zenodo.22798879)


**RAStoERK: A Reference Interaction Atlas for Insights into Signaling and Disease**

![logo](assets/favicon.ico)

* **Authors:** Georg Vucak, Sebastian Didusch, Leandro Cannizzaro, Ana Santiago, Markus Hartl, Jörg Menche, Manuela Baccarini


* **Date:** December 19, 2025


* **Web Access:** rastoerk.univie.ac.at



---

**Core Datasets**

* **Methods Used:** AP-MS and TbID-MS.


* **MPLID23168:** Contains 31 of 35 bait samples along with unincluded extra samples (`data/raw/apms/MPLID23168` and `data/raw/turboid/MPLID23168`).


* **MPLID25011:** Contains 4 remeasured baits processed identically in Spectronaut and post-processing (`data/raw/apms/MPLID25011` and `data/raw/turboid/MPLID25011`).



---

**Directory Structure**

| Folder | Contents |
| --- | --- |
| **data/raw** | Spectronaut output and global analysis files (log2 FC thresholds, bystander list). |
| **data/interim** | Intermediate analysis data structures. |
| **data/processed** | Statistical outputs from limma/DEqMS and Cytoscape `.graphml` files. |
| **data/integration** | IntAct comparison data and downloaded huMAP3.0 complexes. |
| **data/networks** | External network resources. |
| **src/** | Analysis scripts (`src/pipeline` for stats, `src/Figures` for output/tables). |
| **Figures/** | Individual figure panels and source data (`Figure/data`). |

---

**Re-Running the Analysis**

1. Download Spectronaut raw and processed data from rastoerk.univie.ac.at/data.


2. Decompress and extract both folders directly into `data/raw/` and `data/processed/`.


3. For DISGENET gene-disease associations, add your API key to `Figure5_S5/00_pipeline_05_disgenet.R`.
