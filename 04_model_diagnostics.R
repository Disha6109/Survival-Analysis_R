# ==============================================================================
# 04_model_diagnostics.R
#
# PURPOSE
#   Every result in script 03 depends on the PROPORTIONAL HAZARDS assumption:
#   that the hazard ratio for each covariate is constant over time (the two
#   survival curves stay a constant "distance" apart on the hazard scale,
#   rather than crossing or converging). This script checks that assumption
#   instead of just trusting it - this is the step that separates a real
#   biostatistics analysis from a "ran a function, got a p-value" exercise,
#   and it's exactly what a reviewer / interviewer will look for.
#
# CONCEPTUAL BACKGROUND
#   cox.zph() tests the assumption by correlating each covariate's SCALED
#   SCHOENFELD RESIDUALS against time. If a covariate's effect is genuinely
#   constant over time, its residuals should show no trend (a flat line with
#   noise around it). A significant p-value (conventionally < 0.05) or an
#   obvious slope in the diagnostic plot means the PH assumption is violated
#   for that covariate, and the hazard ratio should be interpreted with
#   caution (or the model should be extended, e.g. by stratifying on that
#   covariate or adding a time interaction).
# ==============================================================================

library(survival)
library(survminer)
library(TH.data)

data("GBSG2", package = "TH.data")
surv_obj <- Surv(time = GBSG2$time, event = GBSG2$cens)

cox_model <- coxph(
  surv_obj ~ horTh + age + menostat + tsize + tgrade + pnodes + progrec + estrec,
  data = GBSG2
)

# ---- 1. Formal test of the proportional hazards assumption ------------------
ph_test <- cox.zph(cox_model)
cat("\n--- Proportional hazards test (cox.zph) ---\n")
print(ph_test)

ph_results <- as.data.frame(ph_test$table)
ph_results$covariate <- rownames(ph_results)
rownames(ph_results) <- NULL
ph_results <- ph_results[, c("covariate", "chisq", "df", "p")]
write.csv(ph_results, "outputs/tables/03_ph_assumption_test.csv", row.names = FALSE)

flagged <- ph_results[ph_results$p < 0.05 & ph_results$covariate != "GLOBAL", ]
if (nrow(flagged) > 0) {
  cat("\nCovariates that may VIOLATE the proportional hazards assumption (p < 0.05):\n")
  print(flagged)
} else {
  cat("\nNo covariate significantly violates the proportional hazards assumption.\n")
}

# ---- 2. Diagnostic plots: scaled Schoenfeld residuals vs. time -------------
# A flat line = assumption holds for that covariate. A sloped line = it doesn't.
# plot() on a cox.zph object draws one panel PER COVARIATE - without an
# explicit multi-panel layout, each panel overwrites the last on a single-page
# device, so we set par(mfrow) to a grid big enough for all covariates first.
n_panels <- nrow(ph_test$table) - 1  # exclude the GLOBAL row
n_col <- 3
n_row <- ceiling(n_panels / n_col)

png("outputs/figures/04_schoenfeld_residuals.png",
    width = n_col * 500, height = n_row * 450, res = 130)
par(mfrow = c(n_row, n_col), mar = c(4, 4, 2, 1))
plot(ph_test)
dev.off()
cat("\nSaved outputs/figures/04_schoenfeld_residuals.png\n")

# ---- 3. Influence diagnostics: dfbeta -----------------------------------
# Identifies individual patients who disproportionately influence each
# coefficient estimate - useful for spotting data-entry errors or outlier
# patients before you trust the model.
dfbeta_vals <- residuals(cox_model, type = "dfbeta")
influence_summary <- data.frame(
  covariate      = names(coef(cox_model)),
  max_abs_dfbeta = round(apply(abs(dfbeta_vals), 2, max), 4)
)
write.csv(influence_summary, "outputs/tables/04_influence_diagnostics.csv", row.names = FALSE)
cat("\n--- Max absolute dfbeta per covariate (influence check) ---\n")
print(influence_summary)

# ---- 4. Martingale residuals: check functional form of continuous covariates
# For continuous predictors (age, tsize, pnodes, progrec, estrec), martingale
# residual plots reveal whether the linear term in the model is appropriate,
# or whether the covariate needs a transformation (e.g. log(pnodes + 1)).
null_model <- coxph(surv_obj ~ 1, data = GBSG2)
martingale_res <- residuals(null_model, type = "martingale")

png("outputs/figures/05_martingale_residuals.png", width = 1600, height = 1000, res = 150)
par(mfrow = c(2, 3))
for (var in c("age", "tsize", "pnodes", "progrec", "estrec")) {
  plot(GBSG2[[var]], martingale_res,
       xlab = var, ylab = "Martingale residual",
       main = paste("Martingale residuals vs.", var), pch = 20, col = "#1B9E77")
  lines(lowess(GBSG2[[var]], martingale_res), col = "#E7298A", lwd = 2)
}
dev.off()
cat("Saved outputs/figures/05_martingale_residuals.png\n")

cat("\nDiagnostics complete. Review outputs/tables/03_ph_assumption_test.csv\n")
cat("and outputs/figures/04_schoenfeld_residuals.png before trusting the\n")
cat("hazard ratios reported in script 03.\n")
