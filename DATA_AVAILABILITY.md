# Data Availability

**This repository contains analysis and figure-generation code only. It does not
contain any sequencing data, individual-level (per-donor) data, or clinical
metadata.** The scripts read from local input paths (see each script's header) that
are not distributed here.

## Governance

Several of the datasets analyzed are **controlled-access**. Individual-level data
and donor clinical metadata derived from these sources are **not** redistributed in
this repository and must be obtained directly from the primary source under the
relevant Data Use Agreement (DUA/DUC). Only aggregate / summary results appear in
the manuscript and its supplementary tables.

## Datasets used by this code

Exact accessions, versions, and access procedures are listed in the **Methods**
section of the manuscript. In summary:

| Dataset | Modality | Access | Notes |
|---|---|---|---|
| ROSMAP snRNA-seq (Mathys et al., 2019) | single-nucleus RNA-seq, prefrontal cortex | **Controlled** (RADC / AD Knowledge Portal, DUA required) | Discovery cohort |
| SEA-AD (Allen Institute for Brain Science), middle temporal gyrus | single-nucleus RNA-seq | Open (SEA-AD data portal) | Independent replication cohort |
| Human AD brain bulk RNA-seq (Nativio et al.) | bulk RNA-seq | Public repository (see Methods) | Bulk UPR / GO analysis |
| Additional comparison cohorts referenced in reviewer responses | single-nucleus RNA-seq | See Methods | Used for cross-cohort comparison only |

> Dataset names above identify the primary studies. Do **not** treat this table as an
> access grant — obtain each dataset from its own source under its own terms. For
> controlled data, an approved Data Use Agreement is required before download.

## Reproducing the analysis

1. Obtain the datasets above from their primary sources under the applicable terms.
2. Update the input paths at the top of each script to point to your local copies.
3. Install dependencies (`requirements.txt`; GSVA step needs R/Bioconductor).
4. Run the scripts (see `README.md` for the script → figure map).

No controlled or individual-level data should be committed back into this repository.
