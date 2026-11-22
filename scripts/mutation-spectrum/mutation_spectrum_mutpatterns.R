#!/usr/bin/env Rscript
# =============================================================================
# Mutation spectrum analysis using MutationalPatterns (R/Bioconductor)
# Classifies confirmed TP SNVs into 7-bin spectrum:
#   C>A, C>G, C>T@CpG, C>T_nonCpG, T>A, T>C, T>G
#
# Then compares each somatic (mosaic) spectrum with a pooled germline
# reference: Monte Carlo chi-square goodness-of-fit test, cosine similarity,
# Jensen-Shannon divergence, and pairwise homogeneity tests (BH-adjusted).
#
# Input:  master CSV of confirmed true-positive variants (pre-filtered)
#         + published RAD51C mosaic / germline variant lists
# Output: mutpatterns_output/
#         mutpatterns_output/somatic_germline/
#
# Requirements:
#   BiocManager::install(c("MutationalPatterns","BSgenome.Hsapiens.UCSC.hg38"))
#   install.packages(c("tidyverse","ggplot2"))
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(MutationalPatterns)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(GenomicRanges)
  library(ggplot2)
})

# ---------- Configuration ----------
# Edit these to match your environment. Relative paths resolve against the
# directory this script is launched from; absolute paths also work.
input_csv  <- "all_trios_master.csv"   # master table of manually reviewed variants
output_dir <- "mutpatterns_output"
dir.create(output_dir, showWarnings=FALSE, recursive=TRUE)

# Published RAD51C case: mosaic and het (germline) SNV lists
rad51c_mosaic_csv <- "mosaic-vars-rad51c_truelist.csv"
rad51c_het_csv    <- "het-vars-rad51c_truelist.csv"
somgerm_dir       <- file.path(output_dir, "somatic_germline")

n_monte_carlo <- 1e5     # Monte Carlo replicates for chisq.test / fisher.test
rng_seed      <- 1       # fixed seed so Monte Carlo p-values are reproducible

ref_genome <- BSgenome.Hsapiens.UCSC.hg38

# Trio → sample label mapping (all 5 trios)
trio_map <- tibble(
  trio        = c("brca1","brca2","fanca","fanct","control"),
  sample_name = c("BRCA1","BRCA2","FANCA","FANCT","Control")
)

# 7-bin colors
bin_colors <- c(
  "C>A"        = "#1EBFF0",
  "C>G"        = "#050708",
  "C>T@CpG"    = "#E62725",
  "C>T_nonCpG" = "#FCBD6F",
  "T>A"        = "#CAC9C9",
  "T>C"        = "#A1CF64",
  "T>G"        = "#EDC8C5"
)
bin_order <- names(bin_colors)

# ---------- Step 1: Load & filter master CSV ----------
cat("=== Loading confirmed TP SNVs from master CSV ===\n")

master <- read.csv(input_csv, stringsAsFactors=FALSE)

all_variants <- master |>
  filter(variant_type   == "SNV",
         nchar(REF) == 1, nchar(ALT) == 1,
         REF %in% c("A","C","G","T"),
         ALT %in% c("A","C","G","T")) |>
  left_join(trio_map, by="trio") |>
  dplyr::rename(Chr=CHROM, Ref=REF, Alt=ALT)

cat(sprintf("  %d confirmed TP SNVs across %d samples\n",
            nrow(all_variants), n_distinct(all_variants$sample_name)))
cat("  Per sample:\n")
print(count(all_variants, sample_name))

# ---------- Step 2: Create GRanges ----------
cat("\n=== Building GRanges ===\n")

gr <- GRanges(
  seqnames = all_variants$Chr,
  ranges   = IRanges(start=as.integer(all_variants$POS), width=1),
  ref      = all_variants$Ref,
  alt      = all_variants$Alt,
  sample   = all_variants$sample_name
)
seqlevelsStyle(gr) <- "UCSC"

# ---------- Step 3: Trinucleotide context ----------
cat("\n=== Extracting trinucleotide context ===\n")

trinuc_context <- getSeq(ref_genome, resize(gr, 3, fix="center"))
trinuc_context <- as.character(trinuc_context)
cat(sprintf("  Extracted %d trinucleotide contexts\n", length(trinuc_context)))

# ---------- Step 4: Classify into 7 bins ----------
cat("\n=== Classifying into 7 bins ===\n")

complement <- c(A="T", T="A", C="G", G="C")

classify_7bin <- function(ref, alt, trinuc) {
  if (ref %in% c("A","G")) {
    ref    <- complement[ref]
    alt    <- complement[alt]
    trinuc <- chartr("ACGT","TGCA",trinuc)
    trinuc <- paste0(rev(strsplit(trinuc,"")[[1]]), collapse="")
  }
  change <- paste0(ref, ">", alt)
  if (change == "C>T") {
    if (substr(trinuc, 3, 3) == "G") return("C>T@CpG") else return("C>T_nonCpG")
  }
  return(change)
}

all_variants <- all_variants |>
  mutate(
    trinuc = trinuc_context,
    bin_7  = pmap_chr(list(Ref, Alt, trinuc), classify_7bin),
    bin_7  = factor(bin_7, levels=bin_order)
  )

# ---------- Step 5: Generate 96-channel profile via MutationalPatterns ----------
cat("\n=== Generating MutationalPatterns 96-channel profile ===\n")

vcf_dir <- file.path(output_dir, "vcf_files")
dir.create(vcf_dir, showWarnings=FALSE)

for (sname in trio_map$sample_name) {
  sample_vars <- all_variants |> filter(sample_name == sname)
  vcf_path    <- file.path(vcf_dir, paste0(sname, ".vcf"))
  vcf_lines   <- c(
    "##fileformat=VCFv4.2",
    paste0("##source=confirmed_TP_SNV_", sname),
    paste0("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t", sname)
  )
  data_lines <- sample_vars |>
    mutate(line=paste(Chr, POS, ".", Ref, Alt, ".", "PASS", ".", "GT", "0/1", sep="\t")) |>
    pull(line)
  writeLines(c(vcf_lines, data_lines), vcf_path)
  cat(sprintf("  VCF: %s  (%d SNVs)\n", vcf_path, nrow(sample_vars)))
}

vcf_files         <- file.path(vcf_dir, paste0(trio_map$sample_name, ".vcf"))
names(vcf_files)  <- trio_map$sample_name

# predefined_dbs_mbs=TRUE prevents MutationalPatterns from auto-merging
# adjacent SNVs into DBS variants, preserving all manually reviewed SNVs.
grl_mp <- read_vcfs_as_granges(vcf_files, sample_names=trio_map$sample_name,
                                genome=ref_genome, predefined_dbs_mbs=TRUE)

mut_mat <- mut_matrix(grl_mp, ref_genome)
cat(sprintf("\n  96-channel matrix: %d types x %d samples\n", nrow(mut_mat), ncol(mut_mat)))

write.csv(mut_mat, file.path(output_dir, "sbs96_matrix_mutpatterns.csv"))

pdf(file.path(output_dir, "sbs96_profile_mutpatterns.pdf"), width=14, height=8)
plot_96_profile(mut_mat, condensed=FALSE)
dev.off()
cat("  Saved 96-channel profile: sbs96_profile_mutpatterns.pdf\n")

# ---------- Step 6: Collapse 96 → 7 bins ----------
cat("\n=== Collapsing 96 channels to 7 bins ===\n")

collapse_96_to_7 <- function(mutation_type) {
  parts       <- str_match(mutation_type, "(.+)\\[(.+)>(.+)\\](.+)")
  ref_base    <- parts[,3]
  alt_base    <- parts[,4]
  three_prime <- parts[,5]
  change      <- paste0(ref_base, ">", alt_base)
  case_when(
    change == "C>T" & three_prime == "G" ~ "C>T@CpG",
    change == "C>T"                      ~ "C>T_nonCpG",
    TRUE                                 ~ change
  )
}

seven_bin_mp <- as.data.frame(mut_mat) |>
  rownames_to_column("mutation_type") |>
  mutate(bin_7 = collapse_96_to_7(mutation_type)) |>
  group_by(bin_7) |>
  summarise(across(all_of(trio_map$sample_name), sum), .groups="drop") |>
  mutate(bin_7 = factor(bin_7, levels=bin_order)) |>
  arrange(bin_7)

cat("\n7-bin counts (MutationalPatterns 96-collapsed):\n")
print(as.data.frame(seven_bin_mp))
write_csv(seven_bin_mp, file.path(output_dir, "seven_bin_spectrum_mutpatterns.csv"))

# ---------- Step 7: Manual 7-bin counts ----------
cat("\n=== Manual 7-bin counts (direct from trinuc classification) ===\n")

seven_bin_manual <- all_variants |>
  count(sample_name, bin_7) |>
  pivot_wider(names_from=sample_name, values_from=n, values_fill=0) |>
  arrange(bin_7)

cat("\n7-bin counts (manual):\n")
print(as.data.frame(seven_bin_manual))
write_csv(seven_bin_manual, file.path(output_dir, "seven_bin_spectrum_manual.csv"))

# ---------- Step 8: Plots ----------
cat("\n=== Generating plots ===\n")

plot_7bin <- function(data_long, title, filename) {
  p <- ggplot(data_long, aes(x=bin_7, y=percentage, fill=bin_7)) +
    geom_col(colour="white", linewidth=0.3) +
    geom_text(aes(label=count), vjust=-0.5, size=2.8) +
    facet_wrap(~sample, nrow=1) +
    scale_fill_manual(values=bin_colors) +
    scale_y_continuous(expand=expansion(mult=c(0,0.15))) +
    labs(title=title, x=NULL, y="Percentage (%)") +
    theme_minimal(base_size=11) +
    theme(axis.text.x=element_text(angle=45, hjust=1, size=8),
          legend.position="none",
          strip.text=element_text(face="bold", size=11),
          plot.title=element_text(face="bold", hjust=0.5, size=13),
          panel.grid.major.x=element_blank())
  ggsave(filename, p, width=3.5*n_distinct(data_long$sample), height=5, dpi=150)
  cat(sprintf("  Saved: %s\n", filename))
}

# 96-collapsed 7-bin plot
mp_long <- seven_bin_mp |>
  pivot_longer(cols=all_of(trio_map$sample_name), names_to="sample", values_to="count") |>
  group_by(sample) |>
  mutate(total=sum(count), percentage=count/total*100) |>
  ungroup() |>
  mutate(sample=factor(sample, levels=trio_map$sample_name))

plot_7bin(mp_long,
          "7-Bin Mutation Spectrum — Confirmed TP SNVs (MutationalPatterns)",
          file.path(output_dir, "spectrum_7bin_mutpatterns.pdf"))

# Manual classification plot
manual_long <- seven_bin_manual |>
  pivot_longer(cols=all_of(trio_map$sample_name), names_to="sample", values_to="count") |>
  group_by(sample) |>
  mutate(total=sum(count), percentage=count/total*100) |>
  ungroup() |>
  mutate(sample=factor(sample, levels=trio_map$sample_name))

plot_7bin(manual_long,
          "7-Bin Mutation Spectrum — Confirmed TP SNVs (Manual Classification)",
          file.path(output_dir, "spectrum_7bin_manual.pdf"))

# ---------- Step 9: 96-channel cosine similarity between samples ----------
cat("\n=== Cosine similarity between samples (96-channel) ===\n")

cos_sim_mat <- cos_sim_matrix(mut_mat, mut_mat)
cat("Cosine similarity matrix:\n")
print(round(cos_sim_mat, 3))
write.csv(round(cos_sim_mat, 4), file.path(output_dir, "cosine_similarity_96.csv"))

pdf(file.path(output_dir, "cosine_similarity_heatmap.pdf"), width=6, height=5)
plot_cosine_heatmap(cos_sim_mat, cluster_rows=FALSE, cluster_cols=FALSE)
dev.off()
cat("  Saved: cosine_similarity_heatmap.pdf\n")

# ---------- Step 10: Comparison: MutationalPatterns vs manual ----------
cat("\n=== Comparison: MutationalPatterns vs manual ===\n")

comparison <- seven_bin_mp |>
  pivot_longer(cols=all_of(trio_map$sample_name), names_to="sample", values_to="MP_count") |>
  left_join(
    seven_bin_manual |>
      pivot_longer(cols=all_of(trio_map$sample_name), names_to="sample", values_to="Manual_count"),
    by=c("bin_7","sample")) |>
  mutate(match=MP_count==Manual_count, diff=MP_count-Manual_count)

cat("\nComparison table:\n")
print(as.data.frame(comparison))
write_csv(comparison, file.path(output_dir, "comparison_mp_vs_manual.csv"))

mismatches <- filter(comparison, !match)
if (nrow(mismatches) == 0) {
  cat("\n✓ All counts match between MutationalPatterns and manual classification.\n")
} else {
  cat(sprintf("\n⚠ %d mismatches found:\n", nrow(mismatches)))
  print(as.data.frame(mismatches))
}

# ---------- Step 11: Somatic (mosaic) and germline (het) sets ----------
# Somatic: one sample per individual (RAD51C case + 5 trios, mosaic SNVs).
# Germline reference: het SNVs from the RAD51C case and all 5 trios, pooled.
# Only single-base A/C/G/T SNVs whose REF matches hg38 are kept.
cat("\n=== Step 11: Somatic and germline sets ===\n")
dir.create(somgerm_dir, showWarnings=FALSE, recursive=TRUE)

ref_check_log <- list()

# data.frame -> GRanges, single-base SNVs only, drop sites with REF != hg38
df_to_gr <- function(df, chr_col, pos_col, ref_col, alt_col, label="") {
  bases <- c("A","C","G","T")
  df <- df[nchar(df[[ref_col]])==1 & nchar(df[[alt_col]])==1 &
             df[[ref_col]] %in% bases & df[[alt_col]] %in% bases, , drop=FALSE]
  chr <- as.character(df[[chr_col]])
  if (!any(grepl("^chr", chr))) chr <- paste0("chr", chr)

  gr <- GRanges(
    seqnames = chr,
    ranges   = IRanges(start=as.integer(df[[pos_col]]), width=1),
    strand   = "*",
    REF      = df[[ref_col]],
    ALT      = df[[alt_col]]
  )
  ref_base <- as.character(getSeq(ref_genome, gr))
  ok <- ref_base == as.character(mcols(gr)$REF)
  ref_check_log[[label]] <<- c(snv_in=length(gr), dropped_ref_mismatch=sum(!ok),
                               kept=sum(ok))
  if (any(!ok))
    message("  ", label, ": ", sum(!ok), " site(s) dropped — REF != hg38.")
  gr[ok]
}

# -- 11a. Master CSV, split by consensus_class --------------------------------
tp_snv <- master |> filter(variant_type == "SNV")
master_mosaic <- tp_snv |> filter(consensus_class == "mosaic")
master_het    <- tp_snv |> filter(consensus_class == "het")
cat(sprintf("  Master: %d mosaic + %d het TP SNVs\n",
            nrow(master_mosaic), nrow(master_het)))

# -- 11b. RAD51C published case ----------------------------------------------
# fileEncoding handles the UTF-8 byte-order mark on the 'Identifier' header
rad51c_mosaic_raw <- read.csv(rad51c_mosaic_csv, stringsAsFactors=FALSE,
                              fileEncoding="UTF-8-BOM")

rad51c_het_raw <- read.csv(rad51c_het_csv, stringsAsFactors=FALSE,
                           fileEncoding="UTF-8-BOM")
cat(sprintf("  RAD51C: %d mosaic + %d het SNVs\n",
            nrow(rad51c_mosaic_raw), nrow(rad51c_het_raw)))

# -- 11c. Somatic samples (RAD51C + 5 trios) ---------------------------------
grl_somatic <- c(
  list(RAD51C = df_to_gr(rad51c_mosaic_raw, "Chr","POS","Ref","Alt",
                         label="RAD51C-mosaic")),
  setNames(lapply(trio_map$trio, function(t)
    df_to_gr(master_mosaic[master_mosaic$trio == t, ],
             "CHROM","POS","REF","ALT", label=paste0(t, "-mosaic"))),
    trio_map$sample_name)
)

# -- 11d. Germline reference (RAD51C het + 5 trios het, pooled) --------------
grl_germline <- c(
  list(RAD51C = df_to_gr(rad51c_het_raw, "Chr","POS","Ref","Alt",
                         label="RAD51C-het")),
  setNames(lapply(trio_map$trio, function(t)
    df_to_gr(master_het[master_het$trio == t, ],
             "CHROM","POS","REF","ALT", label=paste0(t, "-het"))),
    trio_map$sample_name)
)
gr_germline_pooled <- do.call(c, unname(grl_germline))

germline_composition <- tibble(
  source = c(names(grl_germline), "TOTAL (Germline reference)"),
  n_snv  = c(sapply(grl_germline, length), length(gr_germline_pooled))
)
cat("\n  Germline reference composition:\n")
print(as.data.frame(germline_composition))
cat(sprintf("  -> Germline reference n = %d\n",
            length(gr_germline_pooled)))
write_csv(germline_composition,
          file.path(somgerm_dir, "germline_reference_composition.csv"))

ref_check_df <- bind_rows(lapply(names(ref_check_log), function(k)
  tibble(set=k, !!!as.list(ref_check_log[[k]]))))
write_csv(ref_check_df, file.path(somgerm_dir, "ref_check_log.csv"))

# ---------- Step 12: 7-bin counts (somatic samples + pooled germline) ----------
cat("\n=== Step 12: 7-bin counts (somatic vs germline) ===\n")

count_7bins <- function(gr) {
  if (length(gr) == 0) return(table(factor(character(0), levels=bin_order)))
  ref    <- as.character(mcols(gr)$REF)
  alt    <- as.character(mcols(gr)$ALT)
  trinuc <- as.character(getSeq(ref_genome, resize(gr, 3, fix="center")))
  bins   <- mapply(classify_7bin, ref, alt, trinuc, USE.NAMES=FALSE)
  table(factor(bins, levels=bin_order))
}

counts7 <- sapply(c(grl_somatic, list(Germline=gr_germline_pooled)), count_7bins)
rownames(counts7) <- bin_order
cat("\n7-bin counts:\n"); print(counts7)
write.csv(counts7, file.path(somgerm_dir, "seven_bin_counts.csv"))

# ---------- Step 13: Goodness-of-fit vs pooled germline ----------
# H0: somatic 7-bin counts ~ Multinomial(n, germline frequencies).
# Pearson chi-square with Monte Carlo p-value (B = n_monte_carlo; smallest
# attainable p ~ 1/(B+1)); BH adjustment across samples.
# Effect sizes: cosine similarity and Jensen-Shannon divergence (nats).
cat("\n=== Step 13: 7-bin GOF vs Germline (Monte Carlo chi-square) ===\n")

p_exp <- counts7[, "Germline"] / sum(counts7[, "Germline"])

cosine_fn <- function(p, q) sum(p*q) / sqrt(sum(p^2) * sum(q^2))
jsd_fn <- function(p, q, eps=1e-12) {
  p <- p/sum(p); q <- q/sum(q); p[p==0] <- eps; q[q==0] <- eps
  m <- 0.5*(p+q)
  0.5*sum(p*log(p/m)) + 0.5*sum(q*log(q/m))           # nats
}

somatic_samples <- setdiff(colnames(counts7), "Germline")
freq7 <- sweep(counts7, 2, colSums(counts7), "/")

set.seed(rng_seed)
gof_res <- bind_rows(lapply(somatic_samples, function(s) {
  obs <- counts7[, s]
  pr  <- freq7[, s]
  pv  <- suppressWarnings(
    chisq.test(obs, p=p_exp, rescale.p=TRUE,
               simulate.p.value=TRUE, B=n_monte_carlo)$p.value)
  tibble(sample          = s,
         n               = sum(obs),
         p_value         = signif(pv, 4),
         cosine_to_germ  = round(cosine_fn(pr, p_exp), 3),
         JSD_to_germ     = round(jsd_fn(pr, p_exp), 4))
})) |>
  mutate(p_adj_BH = signif(p.adjust(p_value, method="BH"), 4)) |>
  relocate(p_adj_BH, .after=p_value)

cat("\nGOF results (vs pooled Germline reference):\n")
print(as.data.frame(arrange(gof_res, p_value)))
write_csv(gof_res, file.path(somgerm_dir, "gof_vs_germline.csv"))

# ---------- Step 14: Pairwise 7-bin homogeneity between somatic samples ----------
# Fisher-Freeman-Halton test on each 2 x 7 table (Monte Carlo p-value),
# BH adjustment across all pairs.
cat("\n=== Step 14: Pairwise 7-bin homogeneity ===\n")

counts7_som <- counts7[, somatic_samples, drop=FALSE]

pairwise_pvals <- function(counts7, B=n_monte_carlo) {
  S <- colnames(counts7); m <- length(S)
  P <- matrix(NA_real_, m, m, dimnames=list(S,S))
  for (i in 1:(m-1)) for (j in (i+1):m) {
    tab <- rbind(counts7[,i], counts7[,j])
    pv  <- suppressWarnings(
      fisher.test(tab, simulate.p.value=TRUE, B=B)$p.value)
    P[i,j] <- P[j,i] <- pv
  }
  diag(P) <- 1; P
}

set.seed(rng_seed)
P <- pairwise_pvals(counts7_som)
upper <- P[upper.tri(P)]
qvec  <- p.adjust(upper, method="BH")
Q <- P; Q[upper.tri(Q)] <- qvec; Q[lower.tri(Q)] <- t(Q)[lower.tri(Q)]; diag(Q) <- 1

cat("\nBH-adjusted q-values (7-bin pairwise homogeneity):\n")
print(round(Q, 4))
write.csv(round(Q, 4), file.path(somgerm_dir, "pairwise_homogeneity_qvals.csv"))

cat(sprintf("\n=== Done. All outputs saved to %s/ ===\n", output_dir))
cat("Outputs:\n")
cat("  sbs96_matrix_mutpatterns.csv   — 96-channel count matrix\n")
cat("  sbs96_profile_mutpatterns.pdf  — MutationalPatterns 96-channel bar plot\n")
cat("  seven_bin_spectrum_mutpatterns.csv / manual.csv — 7-bin counts\n")
cat("  spectrum_7bin_mutpatterns.pdf / manual.pdf      — 7-bin bar plots\n")
cat("  cosine_similarity_96.csv / heatmap.pdf          — sample similarity\n")
cat("  vcf_files/                     — per-sample VCF intermediates\n")
cat("  somatic_germline/germline_reference_composition.csv — n per source (239)\n")
cat("  somatic_germline/ref_check_log.csv                  — REF != hg38 drops\n")
cat("  somatic_germline/seven_bin_counts.csv               — 7-bin, somatic + Germline\n")
cat("  somatic_germline/gof_vs_germline.csv                — MC chi-square GOF, cosine, JSD, BH\n")
cat("  somatic_germline/pairwise_homogeneity_qvals.csv     — pairwise BH q-values\n")
