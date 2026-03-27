# =============================================================================
# CRISPR Target Site Prediction in Plant Genomes
# =============================================================================
# Author      : Alvi (Department of Plant Production & Biotechnology, FAST)
# University  : University of Layyah, Pakistan
# Description : This script identifies and scores CRISPR-Cas9 candidate target
#               sites (protospacers) from a user-supplied DNA sequence, applies
#               on-target scoring heuristics, and flags potential off-target
#               risk based on GC content and seed-region composition.
#
# PAM         : NGG (SpCas9 canonical PAM)
# Guide RNA   : 20 nt protospacer immediately upstream of PAM
# References  :
#   Doench et al. (2016) Nat Biotechnol 34:184-191
#   Hsu et al. (2013) Nat Biotechnol 31:827-832
#   Stemmer et al. (2015) PLoS ONE 10:e0143697
# =============================================================================

# ── 0. Dependencies ──────────────────────────────────────────────────────────
required_packages <- c("stringr", "ggplot2", "dplyr", "tidyr", "Biostrings")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg == "Biostrings") {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install("Biostrings", ask = FALSE)
    } else {
      install.packages(pkg, repos = "https://cran.r-project.org")
    }
  }
}

suppressPackageStartupMessages({
  library(stringr)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(Biostrings)
})

cat("=== CRISPR Target Site Prediction Tool ===\n")
cat("Department of Plant Production & Biotechnology, FAST\n")
cat("University of Layyah, Pakistan\n\n")

# ── 1. Example plant sequence (Arabidopsis thaliana FT gene fragment) ─────────
# Replace this string with your FASTA sequence for real analysis.
# Source: TAIR (https://www.arabidopsis.org), AT1G65480
plant_sequence <- paste0(
  "ATGGATCCAAAGAGGAGTAGCAGCAACAACAACAACAAATGGAGCAGCAGCAGCAGCAGCAGCAGCAGCAG",
  "CAGCAGCAGCAGCAGCAGCAGCAACAGCAGCAGCAGCAAGCAGCAATGGCTTCTTCTTCTTCAGCAGCAGC",
  "AGCAGCAGCAATGGACCCAAAGAGAAGTAGTAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAG",
  "GCTTCCTCCACCGCCGCTGCAGCAATGAACCCAAAGAGAAGCAGCAGCAGCAACAAGCAGCAATGGACCCAA",
  "AGAGAAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAAGCAGCAATGGCTTCCTC",
  "CACCGCCGCAGCAGCAATGAACCCAAAGAGAAGCAGCAACAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGC"
)

cat("Input sequence length:", nchar(plant_sequence), "bp\n\n")

# ── 2. Core Functions ─────────────────────────────────────────────────────────

#' Extract all 20-nt protospacers upstream of NGG PAM on both strands
#'
#' @param sequence Character. Raw DNA sequence (A/T/G/C only, upper-case).
#' @return Data frame with columns: protospacer, position, strand, pam.
find_protospacers <- function(sequence) {
  sequence <- toupper(sequence)
  results  <- list()
  guide_len <- 20L
  pam_len   <-  3L
  window    <- guide_len + pam_len          # 23 nt total

  # Forward strand: look for [20nt][NGG]
  fwd_matches <- gregexpr(
    pattern = paste0("[ATGC]{", guide_len, "}[ATGC]GG"),
    text    = sequence,
    perl    = TRUE
  )[[1]]

  if (fwd_matches[1] != -1) {
    for (pos in fwd_matches) {
      site <- substr(sequence, pos, pos + window - 1)
      proto <- substr(site, 1, guide_len)
      pam   <- substr(site, guide_len + 1, window)
      results[[length(results) + 1]] <- data.frame(
        protospacer = proto,
        position    = pos,
        strand      = "+",
        pam         = pam,
        stringsAsFactors = FALSE
      )
    }
  }

  # Reverse complement strand
  rc_sequence <- as.character(
    reverseComplement(DNAString(sequence))
  )
  seq_len <- nchar(sequence)

  rc_matches <- gregexpr(
    pattern = paste0("[ATGC]{", guide_len, "}[ATGC]GG"),
    text    = rc_sequence,
    perl    = TRUE
  )[[1]]

  if (rc_matches[1] != -1) {
    for (pos in rc_matches) {
      site  <- substr(rc_sequence, pos, pos + window - 1)
      proto <- substr(site, 1, guide_len)
      pam   <- substr(site, guide_len + 1, window)
      # Convert RC position back to forward-strand coordinates
      fwd_pos <- seq_len - (pos + window - 2)
      results[[length(results) + 1]] <- data.frame(
        protospacer = proto,
        position    = fwd_pos,
        strand      = "-",
        pam         = pam,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(results) == 0) {
    cat("No NGG PAM sites detected in the provided sequence.\n")
    return(NULL)
  }

  do.call(rbind, results)
}

#' Calculate GC content of a nucleotide sequence
#'
#' @param seq Character. Nucleotide sequence.
#' @return Numeric. GC fraction (0-1).
gc_content <- function(seq) {
  seq  <- toupper(seq)
  gc   <- str_count(seq, "[GC]")
  gc / nchar(seq)
}

#' Heuristic on-target score (simplified Doench 2016 proxy)
#'
#' This is NOT the full Rule Set 2 model. It approximates key determinants:
#'   - GC content optimum (40-70 %)
#'   - Absence of poly-T (TTTT) which terminates U6 transcription
#'   - Seed region (positions 1-12 from PAM) GC balance
#'
#' @param protospacer Character. 20-nt guide sequence.
#' @return Numeric. Score 0-100 (higher = more favourable).
on_target_score <- function(protospacer) {
  gc     <- gc_content(protospacer)
  seed   <- substr(protospacer, 9, 20)       # positions 9-20 (seed region)
  gc_seed <- gc_content(seed)

  # Penalty for poly-T (terminates Pol III transcription)
  poly_t_penalty <- ifelse(grepl("TTTT", protospacer), 20, 0)

  # GC optimum window: 40-70 % gets full marks, outside loses proportionally
  gc_score <- if (gc >= 0.4 & gc <= 0.7) {
    100
  } else if (gc < 0.4) {
    100 * (gc / 0.4)
  } else {
    100 * ((1 - gc) / 0.3)
  }

  # Seed-region balance bonus (40-60 % GC in seed)
  seed_bonus <- if (gc_seed >= 0.4 & gc_seed <= 0.6) 10 else 0

  score <- max(0, gc_score + seed_bonus - poly_t_penalty)
  round(min(score, 100), 2)
}

#' Flag potential off-target risk based on seed-region composition
#'
#' Guides with high GC in seed region tolerate more mismatches (Hsu et al. 2013).
#'
#' @param protospacer Character. 20-nt guide sequence.
#' @return Character. Risk category: "Low", "Moderate", or "High".
off_target_risk <- function(protospacer) {
  seed    <- substr(protospacer, 9, 20)
  gc_seed <- gc_content(seed)
  if (gc_seed > 0.65) "High" else if (gc_seed >= 0.45) "Moderate" else "Low"
}

# ── 3. Run Prediction ─────────────────────────────────────────────────────────
cat("Scanning sequence for NGG PAM sites...\n")
sites <- find_protospacers(plant_sequence)

if (!is.null(sites)) {
  sites <- sites %>%
    mutate(
      gc_content     = round(sapply(protospacer, gc_content) * 100, 1),
      on_target_score = sapply(protospacer, on_target_score),
      off_target_risk = sapply(protospacer, off_target_risk),
      guide_id        = paste0("gRNA_", seq_len(n()))
    ) %>%
    arrange(desc(on_target_score))

  cat("Total candidate sites identified:", nrow(sites), "\n\n")

  # Top 10 candidates
  top_sites <- head(sites, 10)

  cat("── Top 10 Candidate gRNA Sites ──────────────────────────────────────\n")
  print(
    top_sites %>%
      select(guide_id, protospacer, position, strand, pam,
             gc_content, on_target_score, off_target_risk),
    row.names = FALSE
  )

  # ── 4. Export Results ───────────────────────────────────────────────────────
  output_dir <- "results"
  if (!dir.exists(output_dir)) dir.create(output_dir)

  write.csv(
    sites,
    file      = file.path(output_dir, "all_crispr_sites.csv"),
    row.names = FALSE
  )

  write.csv(
    top_sites,
    file      = file.path(output_dir, "top10_crispr_sites.csv"),
    row.names = FALSE
  )

  cat("\nResults exported to:", output_dir, "\n\n")

  # ── 5. Visualisations ───────────────────────────────────────────────────────

  # 5a. On-target score distribution
  p1 <- ggplot(sites, aes(x = on_target_score, fill = off_target_risk)) +
    geom_histogram(bins = 20, colour = "white", linewidth = 0.3) +
    scale_fill_manual(
      values = c("Low" = "#2ecc71", "Moderate" = "#f39c12", "High" = "#e74c3c"),
      name   = "Off-target Risk"
    ) +
    labs(
      title    = "On-target Score Distribution of Candidate gRNA Sites",
      subtitle = "Plant CRISPR-Cas9 | NGG PAM | SpCas9",
      x        = "On-target Score (0-100)",
      y        = "Number of Sites",
      caption  = "Scoring heuristic adapted from Doench et al. (2016) Nat Biotechnol"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title    = element_text(face = "bold"),
      plot.subtitle = element_text(colour = "grey50"),
      legend.position = "top"
    )

  ggsave(
    file.path(output_dir, "score_distribution.png"),
    plot   = p1,
    width  = 8,
    height = 5,
    dpi    = 300
  )

  # 5b. GC content vs on-target score scatter
  p2 <- ggplot(sites, aes(x = gc_content, y = on_target_score,
                           colour = off_target_risk, shape = strand)) +
    geom_point(size = 3, alpha = 0.8) +
    geom_vline(xintercept = c(40, 70), linetype = "dashed",
               colour = "grey60", linewidth = 0.5) +
    annotate("rect", xmin = 40, xmax = 70, ymin = -Inf, ymax = Inf,
             alpha = 0.05, fill = "#27ae60") +
    annotate("text", x = 55, y = 5, label = "Optimal GC window",
             colour = "#27ae60", size = 3.5, fontface = "italic") +
    scale_colour_manual(
      values = c("Low" = "#2ecc71", "Moderate" = "#f39c12", "High" = "#e74c3c"),
      name   = "Off-target Risk"
    ) +
    labs(
      title    = "GC Content vs. On-target Score",
      subtitle = "Dashed lines mark the 40-70 % optimal GC window",
      x        = "GC Content (%)",
      y        = "On-target Score",
      shape    = "Strand",
      caption  = "Hsu et al. (2013) Nat Biotechnol; Doench et al. (2016) Nat Biotechnol"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title    = element_text(face = "bold"),
      plot.subtitle = element_text(colour = "grey50"),
      legend.position = "right"
    )

  ggsave(
    file.path(output_dir, "gc_vs_score.png"),
    plot   = p2,
    width  = 8,
    height = 5,
    dpi    = 300
  )

  # 5c. Top 10 bar chart
  p3 <- ggplot(top_sites,
               aes(x = reorder(guide_id, on_target_score),
                   y = on_target_score,
                   fill = off_target_risk)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = paste0(on_target_score, "%")),
              hjust = -0.1, size = 3.5) +
    coord_flip(ylim = c(0, 115)) +
    scale_fill_manual(
      values = c("Low" = "#2ecc71", "Moderate" = "#f39c12", "High" = "#e74c3c"),
      name   = "Off-target Risk"
    ) +
    labs(
      title    = "Top 10 Candidate gRNA Sites by On-target Score",
      subtitle = "CRISPR-Cas9 | Plant Genome | NGG PAM",
      x        = NULL,
      y        = "On-target Score",
      caption  = "Dept. of Plant Production & Biotechnology, FAST, University of Layyah"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title    = element_text(face = "bold"),
      plot.subtitle = element_text(colour = "grey50"),
      legend.position = "bottom"
    )

  ggsave(
    file.path(output_dir, "top10_guides.png"),
    plot   = p3,
    width  = 8,
    height = 5,
    dpi    = 300
  )

  cat("Plots saved to:", output_dir, "\n")
  cat("\nDone. Review 'results/' for CSV outputs and PNG figures.\n")
}

