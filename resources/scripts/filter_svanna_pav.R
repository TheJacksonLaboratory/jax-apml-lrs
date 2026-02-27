#!/usr/bin/env Rscript
# filter_svanna_pav.R
# Filters SvAnna-prioritized SV candidates from the PAV SV caller.
#
# Filters applied:
#   - Removes variants with any failed filters
#   - Retains variants with pathogenicity score (psv) >= 4
#
# Usage:
#   Rscript filter_svanna_pav.R <sampleID>_hifiasm_<refID>_pav_svanna.csv
#
# Arguments:
#   <sampleID>_hifiasm_<refID>_pav_svanna.csv  -- SvAnna CSV output file for PAV variants
#
# Output:
#   SvAnna_PAV_Candidates.csv  -- filtered and ranked PAV SV candidates

library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("One argument required: path to the SvAnna CSV output file.\n  Usage: Rscript filter_svanna_pav.R <sampleID>_hifiasm_<refID>_pav_svanna.csv")
}

pav <- read.csv(args[1])

pav <- pav[pav$failed_filters == "", ]
pav <- pav[pav$psv >= 4, ]
pav <- pav %>% mutate(caller = 'PAV')

out_file <- "SvAnna_PAV_Candidates.csv"
write.table(pav, file = out_file, sep = ',', quote = FALSE, row.names = FALSE, col.names = TRUE)
message("Written: ", out_file)