# ADBrainSeq — UPR / ER-stress transcriptional program in Alzheimer's disease

Analysis and figure-generation code for the study of the unfolded protein response
(UPR) / ER-stress transcriptional program across Braak stage in human Alzheimer's
disease (AD) brain, spanning bulk RNA-seq and single-nucleus RNA-seq (snRNA-seq),
with an independent-cohort replication.

> **Status:** manuscript under peer review (Acta Neuropathologica Communications).
> This repository accompanies the paper for code transparency and reproducibility.

---

## ⚠️ Data policy (read first)

This repository ships **code only**. It contains **no sequencing data, no
individual-level (per-donor) data, and no clinical metadata.** Some source datasets
are **controlled-access**; obtain them from their primary sources under the relevant
Data Use Agreement. See [`DATA_AVAILABILITY.md`](DATA_AVAILABILITY.md). A `.gitignore`
safety net blocks common data file types from being committed by accident.

---

## Repository layout

```
code/            analysis + figure-generation scripts (Python; one R script for GSVA)
code/archive/    earlier script iterations, kept for provenance
requirements.txt Python dependencies (R/Bioconductor GSVA is separate)
DATA_AVAILABILITY.md  where to obtain the datasets; governance notes
```

## Script → analysis map

Scripts are grouped by the analysis they produce. Prefixes `r1q*` / `r2q*` denote the
reviewer question each analysis addresses during revision; `build_figure*` scripts
assemble final figure PDFs.

| Prefix | Analysis |
|---|---|
| `r1q1_*` | Bulk RNA-seq UPR gene-set GO / GSEA enrichment and conserved-gene tables |
| `r1q3_*` | snRNA-seq cell-population structure across Braak stage (Figure 3) |
| `r1q5_*` | Cross-modality reproducibility: bulk vs snRNA overlap, heatmaps, Venn |
| `r1q6_*` | Braak-stage correlation of UPR genes |
| `r2q1_*` | UPR sensor transcripts (EIF2AK3/PERK, ERN1/IRE1, ATF6) across Braak, donor-level |
| `r2q2_*` | Donor-level module scores; GSVA / ssGSEA enrichment |
| `r2q4_*` | Covariate / nuclei-count summary tables |
| `r2q6_*` | Cell-type specificity and matched RNA controls |
| `r2q9_*` | Independent-cohort replication (SEA-AD; additional comparison cohorts) |
| `build_figure*`, `build_supp_*` | Final main / supplementary figure assembly (vector PDF) |

## Environment

- Python ≥ 3.10 with the packages in [`requirements.txt`](requirements.txt)
  (`numpy`, `scipy`, `matplotlib`, `openpyxl`, `pypdf`, `python-docx`, `Pillow`,
  `anndata`, `matplotlib-venn`).
- The GSVA enrichment step (`code/r2q2_gsva.R`) requires **R** with the Bioconductor
  package `GSVA`.
- Figure-assembly scripts use the [Poppler](https://poppler.freedesktop.org/) command-line
  tools (`pdftocairo`, `pdftoppm`) for vector PDF handling.

## Reproducing

1. Obtain the datasets listed in [`DATA_AVAILABILITY.md`](DATA_AVAILABILITY.md) from
   their primary sources under the applicable terms.
2. Update the input paths at the top of each script (they point to local copies;
   placeholder roots such as `/data/adbrainseq/...` mark where inputs are expected).
3. Install dependencies, then run the relevant script(s) from the table above.

## Citation

Manuscript under review at *Acta Neuropathologica Communications*. Citation details
will be added on acceptance. Please check back or contact the authors before citing.

## Authors / contact

Goonho Park, PhD — Department of Pathology, Stanford University; VA Palo Alto Health
Care System. Principal Investigator: Jonathan H. Lin.

## License

Licensing is to be determined pending institutional (Stanford University / VA)
review. Until a license is added, all rights are reserved: the code is provided for
review, transparency, and reproducibility. Contact the authors for reuse permissions.
