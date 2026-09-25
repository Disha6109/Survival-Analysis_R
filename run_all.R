# ==============================================================================
# run_all.R
# Runs the full analysis pipeline in order, from the project root.
#   Rscript scripts/run_all.R
# ==============================================================================
message(">> [1/4] Loading and exploring data...")
source("scripts/01_setup_and_data.R")

message("\n>> [2/4] Kaplan-Meier estimation + log-rank test...")
source("scripts/02_kaplan_meier_analysis.R")

message("\n>> [3/4] Cox proportional hazards model...")
source("scripts/03_cox_ph_model.R")

message("\n>> [4/4] Model diagnostics...")
source("scripts/04_model_diagnostics.R")

message("\nAll done. Check outputs/figures/ and outputs/tables/.")
