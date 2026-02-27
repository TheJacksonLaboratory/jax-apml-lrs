#!/usr/bin/env Rscript
# convert_paraphase.R
# Parses a Paraphase JSON output file and extracts copy number (CN) and
# variant information for a set of target genes.
#
# Usage:
#   Rscript convert_paraphase.R <paraphase.json> <targetGenes.txt>
#
# Arguments:
#   paraphase.json   -- JSON output file from Paraphase
#   targetGenes.txt  -- two-column tab-delimited file with columns: GeneSymbol, Chrom
#
# Output:
#   <paraphase.json>_CN.tsv        -- copy number table for target genes
#   <paraphase.json>_variants.tsv  -- haplotype variant table for target genes

library(jsonlite)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Two arguments required. Usage: Rscript convert_paraphase.R <paraphase.json> <targetGenes.txt>")
}

jsonfile       <- args[1]
targetGenes_path <- args[2]

if (!file.exists(jsonfile)) {
  stop("JSON file not found: ", jsonfile)
}
if (!file.exists(targetGenes_path)) {
  stop("Target genes file not found: ", targetGenes_path)
}

targetGenes <- read.table(targetGenes_path, stringsAsFactors = FALSE)
colnames(targetGenes) <- c('GeneSymbol', 'Chrom')
geneChroms <- targetGenes$Chrom
names(geneChroms) <- tolower(targetGenes$GeneSymbol)

jdata <- fromJSON(jsonfile)

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

# Parse a vector of variant strings ("pos_ref_alt") into a data frame
getBedVariants <- function(strvec) {
  if (length(strvec) > 0) {
    s   <- strsplit(strvec, '_')
    pos <- as.numeric(sapply(s, '[[', 1))
    ref <- sapply(s, '[[', 2)
    alt <- sapply(s, '[[', 3)
    return(data.frame(Position = pos, REF = ref, ALT = alt))
  } else {
    return(data.frame(Position = NULL, REF = NULL, ALT = NULL))
  }
}

# Build a variant table across all haplotypes from haplotype_details
getBed_allHapVariants <- function(hapdetails) {
  lvar    <- lapply(hapdetails, function(x) { x$variants })
  poslist <- list()
  for (k in seq_along(lvar)) {
    if (length(lvar[[k]]) > 0) {
      poslist[[names(lvar)[k]]] <- data.frame(
        Hap = names(lvar)[k],
        getBedVariants(lvar[[k]])
      )
    }
  }
  return(bind_rows(poslist))
}

# Extract CN-related fields from a gene entry
makeTable_cn <- function(alist) {
  j <- grep('_cn', names(alist))
  data.frame(Key = names(alist[j]), Value = as.character(alist[j]))
}

# Identify gene entries by all-lowercase names (Paraphase convention)
is_all_lowercase <- function(x) {
  grepl("^[a-z]+$", gsub("\\d", "", x))
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

i_genes <- which(is_all_lowercase(names(jdata)))

# Build and write CN table
cn_list  <- lapply(jdata[i_genes], makeTable_cn)
cn_table <- bind_rows(lapply(seq_along(cn_list), function(i) {
  data.frame(Gene = names(cn_list)[i], CN = cn_list[[i]], stringsAsFactors = FALSE)
}))
cn_out <- paste0(jsonfile, '_CN.tsv')
write.table(cn_table, file = cn_out, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
message("Written: ", cn_out)

# Build and write variant table
var_details <- list()
for (gene in names(jdata)[i_genes]) {
  phas_gene <- jdata[[gene]]
  chr_gene  <- geneChroms[gene]
  if (!is.null(phas_gene$haplotype_details)) {
    pos_df <- getBed_allHapVariants(phas_gene$haplotype_details)
    var_details[[gene]] <- data.frame(
      CHROM = as.character(chr_gene),
      COORD = pos_df$Position,
      GENE  = gene,
      HAP   = pos_df$Hap,
      REF   = pos_df$REF,
      ALT   = pos_df$ALT
    )
  }
}
var_out <- paste0(jsonfile, '_variants.tsv')
write.table(bind_rows(var_details), file = var_out, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
message("Written: ", var_out)