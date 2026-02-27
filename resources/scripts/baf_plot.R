# baf_plot.R
# Generates B-allele frequency (BAF) and copy number (CN) plots for the JAX-APML LRS pipeline.
#
# The B-allele frequency is derived from the phased VCF output of WhatsHap + DeepVariant.
# Only variants passing the PASS filter are included.
# Copy number data is taken from the HiFiCNV bedgraph output.
#
# Usage:
#   Rscript baf_plot.R <SAMPLE_ID> [<REFERENCE_GENOME_ID>]
#
# Arguments:
#   SAMPLE_ID             -- sample identifier (must match output filenames from pipeline)
#   REFERENCE_GENOME_ID   -- genome build identifier (default: hg38)
#
# Input files (expected in working directory):
#   <SAMPLE_ID>_<REF>_pbmm2_deepvariant_Emedgene.vcf  -- phased DeepVariant VCF
#   hificnv.<SAMPLE_ID>.copynum.bedgraph               -- HiFiCNV copy number bedgraph
#
# Output:
#   <SAMPLE_ID>_<REF>_B-Allele_Plot.pdf  -- per-chromosome BAF and CN plots

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("No arguments provided. Usage: Rscript baf_plot.R <SAMPLE_ID> [<REFERENCE_GENOME_ID>]")
}

sID <- args[1]

if (length(args) < 2) {
  genref <- 'hg38'
  message('No genome build provided. Using default: ', genref)
} else {
  genref <- args[2]
}

cat(version$version.string, "\n")
suppressMessages({
  library(CopyNumberPlots)
  library(VariantAnnotation)
})

vcf_file <- paste0('./', sID, '_', genref, '_pbmm2_deepvariant_Emedgene.vcf')
bdg_file <- paste0('./hificnv.', sID, '.copynum.bedgraph')
chroms   <- paste0('chr', c(1:22, 'X', 'Y'))

message("Reading ", vcf_file)
vcf          <- readVcf(vcf_file, genome = genref)
pass_records <- vcf[rowRanges(vcf)$FILTER == "PASS"]
vaf_data     <- geno(pass_records)$VAF
vaf_sums     <- sapply(vaf_data, sum, na.rm = TRUE)

message("Reading ", bdg_file)
bdg_gr <- import(bdg_file, format = 'bedGraph')
vaf_gr <- rowRanges(pass_records)

# Filter to chromosomes of interest and remove NA values
vaf_gr   <- vaf_gr[seqnames(vaf_gr) %in% chroms & !is.na(vaf_sums)]
vaf_sums <- vaf_sums[!is.na(vaf_sums)]
bdg_gr   <- bdg_gr[seqnames(bdg_gr) %in% chroms & !is.na(mcols(bdg_gr)$score)]

# Overlap VAF positions with CN segments
overlaps   <- findOverlaps(vaf_gr, bdg_gr)
results_gr <- GRanges(
  seqnames = seqnames(vaf_gr[queryHits(overlaps)]),
  ranges   = IRanges(
    start = start(vaf_gr[queryHits(overlaps)]),
    end   = end(vaf_gr[queryHits(overlaps)])
  ),
  vaf = vaf_sums[queryHits(overlaps)],
  cn  = mcols(bdg_gr)$score[subjectHits(overlaps)]
)
results_gr$b_allele_CN <- results_gr$vaf * results_gr$cn

# ---------------------------------------------------------------------------
# Per-chromosome BAF + CN karyotype plot
# ---------------------------------------------------------------------------
plot_vafCN <- function(gr, selCHR, r0panel2 = 0.15,
                       colpoint = '#1F77B488',
                       bg1col   = '#fcfcf5',
                       bg2col   = '#fcfcf5',
                       colbcn   = '#23dea9') {

  agr <- subset(gr, seqnames == selCHR)
  bgr <- subset(agr, cn != 2) # non-diploid subset

  kp <- plotKaryotype(plot.type = 2, chromosomes = selCHR, genome = genref)
  title(main = sID)
  kpAddCytobandLabels(kp)
  kpDataBackground(kp, data.panel = 1, col = bg1col)
  kpDataBackground(kp, data.panel = 2, col = bg2col, r0 = r0panel2)

  kpAddLabels(kp, labels = "B-ratio",  data.panel = 1, cex = 0.7)
  kpAddLabels(kp, labels = "Total CN", data.panel = 2, cex = 0.7, col = 'grey')
  kpAddLabels(kp, labels = "B-CN",     data.panel = 2, cex = 0.7, col = colbcn, r1 = 0.6)

  kpAxis(kp, data.panel = 1, side = 2, cex = 0.7)
  kpAxis(kp, data.panel = 2, side = 2, cex = 0.7,
         ymin = 0, ymax = max(agr$cn),
         tick.pos = c(0, 1, 2, max(agr$cn)), r0 = r0panel2)
  kpAddBaseNumbers(kp, tick.len = 10, tick.col = "navy", cex = 0.7, minor.tick.col = "gray")

  # B-allele frequency points
  kpPoints(kp, data = agr, y = agr$vaf, data.panel = 1, col = colpoint, cex = 0.2)

  # Copy number
  kpAbline(kp, h = 2, lty = 2, r0 = r0panel2, data.panel = 2,
           col = '#c0c0c0', ymin = 0, ymax = max(agr$cn))
  kpPoints(kp, data = agr, y = agr$b_allele_CN, data.panel = 2,
           col = colbcn, cex = 0.2,
           ymax = max(agr$b_allele_CN), ymin = min(agr$b_allele_CN), r0 = r0panel2)

  if (length(bgr) > 0) {
    kpPoints(kp, data = bgr, y = bgr$cn, data.panel = 2,
             ymin = 0, ymax = max(agr$cn),
             pch = ifelse(bgr$cn > 2, 2, 6),
             col = ifelse(bgr$cn > 2, '#1ebbe3bb', '#de3623bb'),
             r0 = r0panel2, cex = 0.8)
  }

  return(agr)
}

# ---------------------------------------------------------------------------
# Generate PDF
# ---------------------------------------------------------------------------
out_pdf <- paste0(sID, '_', genref, '_B-Allele_Plot.pdf')
pdf(out_pdf, width = 10, height = 6)
for (ch in chroms) {
  plot_vafCN(gr = results_gr, selCHR = ch)
  message('Plotting complete for ', ch)
}
dev.off()
message(out_pdf, ' was generated.')