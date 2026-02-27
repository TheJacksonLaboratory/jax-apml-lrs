#!/usr/bin/env Rscript
# filter_svanna_pbsv_sniffles2.R
# Merges and filters SvAnna-prioritized SV candidates from pbsv and Sniffles2 callers.
#
# Filters applied:
#   - Removes variants with any failed filters
#   - Retains variants with pathogenicity score (psv) >= 2
#   - Removes BND (breakend) calls from Sniffles2
#   - Sorts by psv score (descending)
#
# Usage:
#   Rscript filter_svanna_pbsv_sniffles2.R <sampleID>_<refID>_pbmm2_pbsv_svanna.csv <sampleID>_<refID>_pbmm2_sniffles2_svanna.csv
#
# Arguments:
#   <sampleID>_<refID>_pbmm2_pbsv_svanna.csv      -- SvAnna CSV output for PBSV variants
#   <sampleID>_<refID>_pbmm2_sniffles2_svanna.csv  -- SvAnna CSV output for Sniffles2 variants
#
# Output:
#   SvAnna_SV_Candidates.csv  -- merged, filtered, and ranked SV candidates

library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Two arguments required: paths to the PBSV and Sniffles2 SvAnna CSV output files.\n  Usage: Rscript filter_svanna_pbsv_sniffles2.R <sampleID>_<refID>_pbmm2_pbsv_svanna.csv <sampleID>_<refID>_pbmm2_sniffles2_svanna.csv")
}

pbsv      <- read.csv(args[1])
sniffles2 <- read.csv(args[2])

# Filter PBSV variants
pbsv <- pbsv[pbsv$failed_filters == "", ]
pbsv <- pbsv[pbsv$psv >= 2, ]
pbsv <- pbsv %>% mutate(caller = 'pbsv')

# Filter Sniffles2 variants
sniffles2 <- sniffles2[sniffles2$failed_filters == "", ]
sniffles2 <- sniffles2[sniffles2$psv >= 2, ]
sniffles2 <- sniffles2[sniffles2$vtype != "BND", ]
sniffles2 <- sniffles2 %>% mutate(caller = 'Sniffles2')

# Merge and rank
sv <- bind_rows(sniffles2, pbsv)
sv <- sv[order(sv$psv, decreasing = TRUE), ]

out_file <- "SvAnna_SV_Candidates.csv"
write.table(sv, file = out_file, sep = ',', quote = FALSE, row.names = FALSE, col.names = TRUE)
message("Written: ", out_file)