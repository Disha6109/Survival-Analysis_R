# ==============================================================================
# 03_cox_ph_model.R
#
# PURPOSE
#   Kaplan-Meier + log-rank (script 02) can only compare survival between
#   groups of a SINGLE categorical variable, and can't adjust for confounders.
#   The Cox proportional hazards model extends this: it estimates how MULTIPLE
#   covariates (treatment, age, tumor grade, lymph node count, receptor
#   status...) simultaneously affect the hazard of recurrence, while leaving
#   the baseline hazard shape unspecified (that's what makes it "semi-
#   parametric" - a major reason it's the default tool in clinical trial
#   biostatistics).
#
# CONCEPTUAL BACKGROUND
#   The hazard h(t) is the instantaneous risk of the event at time t, given
#   survival up to t. The Cox model assumes:
#       h(t | X) = h0(t) * exp(b1*X1 + b2*X2 + ... + bk*Xk)
#   h0(t) is an unspecified baseline hazard (common to everyone), and the
#   exp(...) term scales that baseline up or down depending on a patient's
#   covariates. exp(b_i) is the HAZARD RATIO for a one-unit increase in X_i,
#   holding everything else constant - e.g. HR = 1.5 for a covariate means a
#   57% higher instantaneous risk of recurrence than someone with a hazard
#   ratio of ~0.95 for that same predictor.  A HAZARD RATIO of 1.5 for a
#   binary covariate means that group has 50% higher instantaneous risk of the
#   event than the reference group, at every point in time (that "every point
#   in time" bit is the PROPORTIONAL HAZARDS assumption - checked in script
#   04).
# ==============================================================================

library(survival)
library(survminer)
library(TH.data)

data("GBSG2", package = "TH.data")
surv_obj <- Surv(time = GBSG2$time, event = GBSG2$cens)

# ---- 1. Multivariable Cox model ---------------------------------------------
# We adjust the treatment effect (horTh) for the standard clinical covariates
# reported in the original Schumacher et al. (1994) analysis: age, menopausal
# status, tumor size, tumor grade, number of positive nodes, and hormone
# receptor levels.
cox_model <- coxph(
  surv_obj ~ horTh + age + menostat + tsize + tgrade + pnodes + progrec + estrec,
  data = GBSG2
)

cat("\n--- Cox proportional hazards model summary ---\n")
print(summary(cox_model))

# ---- 2. Extract hazard ratios + 95% CIs into a clean table -----------------
cox_summary   <- summary(cox_model)
hazard_ratios <- data.frame(
  covariate     = rownames(cox_summary$coefficients),
  hazard_ratio  = round(cox_summary$coefficients[, "exp(coef)"], 3),
  ci_lower      = round(cox_summary$conf.int[, "lower .95"], 3),
  ci_upper      = round(cox_summary$conf.int[, "upper .95"], 3),
  p_value       = round(cox_summary$coefficients[, "Pr(>|z|)"], 4)
)
rownames(hazard_ratios) <- NULL

cat("\n--- Hazard ratios (95% CI) ---\n")
print(hazard_ratios)

write.csv(hazard_ratios, "outputs/tables/02_cox_hazard_ratios.csv", row.names = FALSE)

# ---- 3. Forest plot of hazard ratios ----------------------------------------
forest_plot <- ggforest(cox_model, data = GBSG2)
ggsave(
  filename = "outputs/figures/03_cox_forest_plot.png",
  plot     = forest_plot,
  width    = 8, height = 6, dpi = 300
)
cat("\nSaved outputs/figures/03_cox_forest_plot.png\n")

# ---- 4. Concordance index (model discrimination) ----------------------------
# The C-index is survival analysis's version of AUC: the probability that,
# for a random pair of patients, the model correctly ranks who recurs first.
# 0.5 = no better than chance, 1.0 = perfect discrimination.
cat(sprintf(
  "\nConcordance (C-index): %.3f (se = %.3f)\n",
  cox_summary$concordance["C"], cox_summary$concordance["se(C)"]
))

# ---- 5. Interpretation cheat-sheet, saved alongside the results ------------
interpretation <- c(
  "How to read the hazard ratios in 02_cox_hazard_ratios.csv:",
  "- HR = 1   -> covariate has no effect on recurrence risk",
  "- HR > 1   -> covariate INCREASES the hazard (worse prognosis)",
  "- HR < 1   -> covariate DECREASES the hazard (protective/better prognosis)",
  "- If the 95% CI for HR crosses 1, the effect is not statistically significant",
  "",
  "Example: horThyes HR < 1 would mean hormonal therapy patients have a lower",
  "recurrence hazard than the no-therapy reference group, after adjusting for",
  "age, tumor size/grade, nodes, and receptor levels."
)
writeLines(interpretation, "outputs/tables/02_how_to_read_hazard_ratios.txt")
