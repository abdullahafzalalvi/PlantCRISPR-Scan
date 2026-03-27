# CRISPR Target Site Prediction in Plant Genomes

[![R](https://img.shields.io/badge/R-%3E%3D4.1.0-blue?logo=r)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Bioconductor](https://img.shields.io/badge/Bioconductor-Biostrings-orange)](https://bioconductor.org/packages/Biostrings/)
[![Department](https://img.shields.io/badge/Dept.-Plant%20Production%20%26%20Biotechnology-brightgreen)](https://uol.edu.pk/)

> A lightweight R-based pipeline for identifying and scoring CRISPR-Cas9 candidate target sites (protospacers) from plant genomic DNA sequences, with on-target scoring heuristics and off-target risk flagging.

---

## Overview

CRISPR-Cas9 genome editing requires the identification of suitable 20-nucleotide protospacer sequences upstream of a canonical NGG PAM (protospacer adjacent motif). Selecting high-quality guide RNAs (gRNAs) with favourable on-target activity and low off-target potential is a critical bottleneck in plant genome engineering workflows (Doench et al., 2016; Hsu et al., 2013).

This tool was developed as part of coursework and research activities in the **Department of Plant Production & Biotechnology, Faculty of Agricultural Sciences and Technology (FAST), University of Layyah, Pakistan**. It is intended as an educational and exploratory resource — not a replacement for validated tools such as CRISPOR or Cas-OFFinder in production workflows.

---

## Features

- Scans both forward (+) and reverse complement (−) strands for NGG PAM sites
- Extracts 20-nt protospacer sequences for SpCas9
- Computes GC content per guide
- Applies a heuristic on-target score (0–100) based on:
  - GC content optimum (40–70 %)
  - Seed-region GC balance (positions 9–20 relative to PAM)
  - Poly-T penalty (TTTT terminates U6 Pol III transcription)
- Flags off-target risk (Low / Moderate / High) based on seed-region GC
- Exports results to CSV
- Generates publication-ready PNG plots (score distribution, GC vs. score scatter, top-10 bar chart)

---

## Biological Rationale

| Parameter | Basis | Reference |
|---|---|---|
| GC content 40–70 % | Empirically associated with efficient Cas9 cleavage | Doench et al. (2016) |
| Seed region (PAM-proximal 12 nt) | Mismatch tolerance highest in distal region | Hsu et al. (2013) |
| Poly-T avoidance | TTTT terminates RNA Pol III (U6 promoter) transcription | Stemmer et al. (2015) |
| Off-target risk proxy | Seed GC > 65 % correlates with higher mismatch tolerance | Hsu et al. (2013) |

> **Important caveat:** The scoring function here is a simplified heuristic, not the full Doench Rule Set 2 machine-learning model. For rigorous experimental design, users should validate shortlisted guides using [CRISPOR](http://crispor.tefor.net/) or [Cas-OFFinder](http://www.rgenome.net/cas-offinder/).

---

## Requirements

- R >= 4.1.0
- CRAN packages: `stringr`, `ggplot2`, `dplyr`, `tidyr`
- Bioconductor package: `Biostrings`

The script auto-installs all dependencies on first run.

---

## Installation

```bash
# Clone the repository
git clone https://github.com/<your-username>/crispr-plant-target.git
cd crispr-plant-target
```

---

## Usage

1. Open `crispr_target_prediction.R` in RStudio or any R environment.
2. Replace the `plant_sequence` variable (line ~50) with your target DNA sequence (plain string, A/T/G/C only, any length).
3. Run the script:

```r
source("crispr_target_prediction.R")
```

4. Outputs will appear in the `results/` directory:

```
results/
├── all_crispr_sites.csv      # All detected candidate sites
├── top10_crispr_sites.csv    # Top 10 by on-target score
├── score_distribution.png    # Histogram of score distribution
├── gc_vs_score.png           # GC content vs. on-target score scatter
└── top10_guides.png          # Bar chart of top 10 candidates
```

---

## Example Output

### Console summary

```
=== CRISPR Target Site Prediction Tool ===
Department of Plant Production & Biotechnology, FAST
University of Layyah, Pakistan

Input sequence length: 420 bp
Scanning sequence for NGG PAM sites...
Total candidate sites identified: 38

── Top 10 Candidate gRNA Sites ──────────────────────────────────────
 guide_id          protospacer  position strand pam gc_content on_target_score off_target_risk
  gRNA_1  GCAGCAGCAGCAGCAGCAGC       145      +  CGG       65.0           95.5        Moderate
  gRNA_2  ATGGCTTCTTCTTCTTCAGC       201      -  AGG       45.0           90.0             Low
  ...
```

---

## Limitations

1. This pipeline does not perform genome-wide off-target analysis. Use Cas-OFFinder for exhaustive off-target search against full reference genomes.
2. The on-target score is a heuristic proxy, not a machine-learning prediction. Experimental validation (e.g., T7E1 assay, Sanger sequencing) is mandatory before concluding guide activity.
3. The tool currently supports SpCas9 (NGG PAM) only. Other variants (SaCas9: NNGRRT; Cas12a: TTTV) are not implemented.
4. Input sequences should be free of ambiguous IUPAC codes (N, R, Y, etc.) for reliable PAM detection.

---

## References

- Doench, J. G., Fusi, N., Sullender, M., Hegde, M., Vaimberg, E. W., Donovan, K. F., ... & Root, D. E. (2016). Optimized sgRNA design to maximize activity and minimize off-target effects of CRISPR-Cas9. *Nature Biotechnology*, 34(2), 184–191. https://doi.org/10.1038/nbt.3437
- Hsu, P. D., Scott, D. A., Weinstein, J. A., Ran, F. A., Konermann, S., Agarwala, V., ... & Zhang, F. (2013). DNA targeting specificity of RNA-guided Cas9 nucleases. *Nature Biotechnology*, 31(9), 827–832. https://doi.org/10.1038/nbt.2647
- Stemmer, M., Thumberger, T., del Sol Keyer, M., Wittbrodt, J., & Mateo, J. L. (2015). CCTop: an intuitive, flexible and reliable CRISPR/Cas9 target prediction tool. *PLoS ONE*, 10(4), e0124633. https://doi.org/10.1371/journal.pone.0124633
- The Arabidopsis Information Resource (TAIR). https://www.arabidopsis.org

---

## License

This project is released under the [MIT License](LICENSE).

---

## Author

**Abdullah Afzal Alvi**  
Undergraduate Researcher, 6th Semester  
Department of Plant Production & Biotechnology  
Faculty of Agricultural Sciences and Technology (FAST)  
University of Layyah, Pakistan  

*Supervised by: Dr. Zeshan Hassan*

---

## Acknowledgements

Example sequence used for demonstration is derived from a publicly available *Arabidopsis thaliana* gene fragment (TAIR: AT1G65480). No proprietary data are included in this repository.
