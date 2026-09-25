# ==============================================================================
# 01_setup_and_data.R
#
# PURPOSE
#   Load the GBSG2 dataset (German Breast Cancer Study Group, 1994), inspect
#   it, and write a copy to data/gbsg2.csv so the project has a plain-text
#   dataset sitting in the repo (not just an R package object). This is the
#   dataset every later script reads from.
#
# WHY THIS DATASET
#   GBSG2 is a real clinical trial dataset: 686 women with primary node-
#   positive breast cancer, followed for recurrence-free survival and overall
#   survival after adjuvant chemotherapy. It's a standard teaching/benchmark
#   dataset in survival analysis (used in Schumacher et al., 1994), which is
#   exactly why it's a safe, well-documented choice for a portfolio project:
#   anyone reviewing your repo can independently verify your numbers against
#   the published literature.
#
# WHERE IT COMES FROM
#   Ships inside the CRAN package `TH.data` (Torsten Hothorn's teaching data
#   package) - no manual download, no license issues, fully reproducible from
#   a fresh R install.
# ==============================================================================

# ---- 1. Load packages -------------------------------------------------------
# survival : core time-to-event modelling (Surv, survfit, coxph, cox.zph)
# TH.data  : ships the GBSG2 dataset
library(survival)
library(TH.data)

# ---- 2. Load the dataset -----------------------------------------------------
data("GBSG2", package = "TH.data")

# ---- 3. What are we actually looking at? ------------------------------------
# GBSG2 columns (686 rows, one per patient):
#   horTh     : hormonal therapy - "yes"/"no"                (the TREATMENT arm)
#   age       : age in years
#   menostat  : menopausal status - "Pre"/"Post"
#   tsize     : tumor size (mm)
#   tgrade    : tumor grade - I / II / III (histological grade, ordered factor)
#   pnodes    : number of positive lymph nodes
#   progrec   : progesterone receptor (fmol/l)
#   estrec    : estrogen receptor (fmol/l)
#   time      : recurrence-free survival time, in DAYS               <- time
#   cens      : event indicator - 1 = recurrence/death observed,
#                                  0 = censored (patient event-free at last
#                                      follow-up, or lost to follow-up)  <- status

cat("\n--- Dataset dimensions ---\n")
print(dim(GBSG2))

cat("\n--- Structure ---\n")
str(GBSG2)

cat("\n--- First rows ---\n")
print(head(GBSG2))

cat("\n--- Event / censoring balance ---\n")
# How many patients had the event (recurrence) vs were censored?
print(table(Event = GBSG2$cens, Label = ifelse(GBSG2$cens == 1, "event", "censored")))

cat("\n--- Treatment arm sizes ---\n")
print(table(GBSG2$horTh))

cat("\n--- Missingness check ---\n")
print(colSums(is.na(GBSG2)))

# ---- 4. Save a plain-text copy of the dataset into the repo -----------------
# This means anyone cloning the repo has the exact data file used, without
# needing to know it's bundled inside an R package.
write.csv(GBSG2, "data/gbsg2.csv", row.names = FALSE)
cat("\nSaved data/gbsg2.csv (", nrow(GBSG2), "rows x", ncol(GBSG2), "cols )\n")

# ---- 5. Quick descriptive summary table, saved for the README --------------
summary_stats <- data.frame(
  n_patients        = nrow(GBSG2),
  n_events          = sum(GBSG2$cens == 1),
  n_censored        = sum(GBSG2$cens == 0),
  pct_event         = round(100 * mean(GBSG2$cens == 1), 1),
  median_age        = median(GBSG2$age),
  median_followup_d = median(GBSG2$time)
)
write.csv(summary_stats, "outputs/tables/00_dataset_summary.csv", row.names = FALSE)
print(summary_stats)
