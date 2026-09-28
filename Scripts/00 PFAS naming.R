



# 00_PFAS_Naming_Configuration.R

# This initial configuration code is for PFAS compound naming. 
# Column names get cleaned right at the start of the script. 

# Nothing here reads or writes files.

# ---- 1.0 Column-naming convention ------------------------------------------------
# Cleaner shared by the results table and the MDL table (script 03).
# Result columns look like "4:2FTS Results" -> "X4.2FTS_Results".
# The leading "X" (added by make.names) is KEPT on purpose: names such as
# "4.2FTS" are not valid R names and cannot be used on the left side of a
# model formula, but "X4.2FTS" can.
clean_names <- function(x) {
  x |>
    trimws() |>
    gsub(" ", "_", x = _) |>
    gsub("[-:]", ".", x = _) |>
    make.names()
}

# ---- 2.0 PFAS compounds and functional groups ------------------------------------
# Single source of PFAS naming conventions. 
# script 03 will sums these into group columns, script 04 will model them, 
# script 05 colours plots by them.
pfas_groups <- list(
  PFCA         = c("PFBA", "PFPeA", "PFHxA", "PFHpA", "PFOA",
                   "PFNA", "PFDA", "PFUnA", "PFDoA", "PFTrDA", "PFTeDA"),
  PFSA         = c("PFBS", "PFPeS", "PFHxS", "PFHpS", "PFOS",
                   "PFNS", "PFDS", "PFDoS"),
  FTSA         = c("X4.2FTS", "X6.2FTS", "X8.2FTS"),
  Sulfonamides = c("PFOSA", "NMeFOSA", "NEtFOSA", "NMeFOSE", "NEtFOSE"),
  FOSAA        = c("NMeFOSAA", "NEtFOSAA"),
  PFECA        = c("HFPO.DA", "ADONA", "PFMPA", "PFMBA", "NFDHA"),
  PFESA        = c("X9Cl.PF3ONS", "X11Cl.PF3OUdS", "PFEESA"),
  FTCA         = c("X3.3_FTCA", "X5.3_FTCA", "X7.3_FTCA"))

pfas_individual <- unlist(pfas_groups, use.names = FALSE)   # the 40 measured compounds
pfas_compounds  <- c("PFAS40", pfas_individual, names(pfas_groups))  # every response variable

# Family lookup used to colour/group plots in script 05 (built from pfas_groups
# so plots always agree with the group sums that were modeled)
pfas_family_lookup <- dplyr::bind_rows(
  tibble::tibble(compound = c("PFAS40", names(pfas_groups)),
                 family   = "Summary total"),
  tibble::tibble(family   = rep(names(pfas_groups), lengths(pfas_groups)),
                 compound = pfas_individual))

stopifnot(length(pfas_individual) == 40, !anyDuplicated(pfas_compounds))

