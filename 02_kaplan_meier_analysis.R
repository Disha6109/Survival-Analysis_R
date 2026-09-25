# ==============================================================================
# 02_kaplan_meier_analysis.R
#
# PURPOSE
#   Estimate recurrence-free survival with the Kaplan-Meier method, compare
#   the hormonal-therapy arm against the no-therapy arm, and test whether the
#   difference is statistically significant with a log-rank test.
#
# CONCEPTUAL BACKGROUND (read this before the code)
#   Survival data is different from ordinary numeric data because of
#   CENSORING: some patients hadn't had a recurrence by the time the study
#   ended, or dropped out. We don't know their true survival time - only that
#   it's AT LEAST as long as their observed follow-up. Throwing those patients
#   out (or treating them as "cured") would bias the analysis. Kaplan-Meier
#   handles this correctly: at every time a recurrence happens, it computes
#   the probability of surviving that instant given everyone still at risk,
#   and multiplies those conditional probabilities together to get a step
#   function, S(t) = P(surviving past time t). Patients who are censored are
#   removed from the "at risk" set at the point they're censored, but they
#   still contribute information up to that point.
#
#   The log-rank test then asks: are the two arms' KM curves different enough
#   that it's unlikely to be chance? It compares observed vs. expected events
#   in each arm at every event time, summed up into a chi-squared statistic.
# ==============================================================================

library(survival)
library(survminer)   # ggplot2-based survival curve plotting
library(TH.data)

data("GBSG2", package = "TH.data")

# ---- 1. Build the survival object -------------------------------------------
# Surv(time, event) packages the two columns every survival function needs:
#   time = follow-up duration, cens = 1 if the event (recurrence) was observed
surv_obj <- Surv(time = GBSG2$time, event = GBSG2$cens)

# ---- 2. Overall Kaplan-Meier curve (no grouping) ----------------------------
km_overall <- survfit(surv_obj ~ 1, data = GBSG2)
cat("\n--- Overall KM summary (median survival, 95% CI) ---\n")
print(km_overall)

# ---- 3. Kaplan-Meier curves BY TREATMENT ARM (hormonal therapy) ------------
km_by_arm <- survfit(surv_obj ~ horTh, data = GBSG2)
cat("\n--- KM summary by hormonal therapy arm ---\n")
print(km_by_arm)

# ---- 4. Log-rank test: is the difference between arms significant? --------
logrank_test <- survdiff(surv_obj ~ horTh, data = GBSG2)
cat("\n--- Log-rank test: hormonal therapy vs. no therapy ---\n")
print(logrank_test)

# Pull out a clean p-value for reporting
logrank_p <- 1 - pchisq(logrank_test$chisq, length(logrank_test$n) - 1)
cat(sprintf("\nLog-rank p-value: %.4f\n", logrank_p))

# ---- 5. Plot: Kaplan-Meier curves with risk table and p-value --------------
km_plot <- ggsurvplot(
  km_by_arm,
  data          = GBSG2,
  pval          = TRUE,               # overlays the log-rank p-value
  conf.int      = TRUE,                # 95% confidence bands
  risk.table    = TRUE,                # number-at-risk table beneath the plot
  risk.table.col = "strata",
  legend.labs   = c("No hormonal therapy", "Hormonal therapy"),
  legend.title  = "Treatment arm",
  xlab          = "Time (days)",
  ylab          = "Recurrence-free survival probability",
  title         = "Kaplan-Meier Curves: GBSG2 Breast Cancer Trial",
  palette       = c("#E7298A", "#1B9E77"),
  ggtheme       = theme_minimal()
)

# ggsurvplot() returns a special "ggsurvplot" object (a list containing the
# main plot + risk table as separate ggplot grobs) - NOT a single ggplot, so
# ggsave(plot = ...) can't render it directly. The reliable way to save both
# panels together is to open a graphics device, print() the object into it,
# then close the device.
png("outputs/figures/01_kaplan_meier_curves.png", width = 8, height = 7,
    units = "in", res = 300)
print(km_plot)
dev.off()

cat("\nSaved outputs/figures/01_kaplan_meier_curves.png\n")

# ---- 6. Also stratify by tumor grade, out of clinical interest -------------
# (secondary, exploratory - not the primary comparison of the project, but
#  shows the KM method generalizes to any categorical grouping variable)
km_by_grade <- survfit(surv_obj ~ tgrade, data = GBSG2)
grade_plot <- ggsurvplot(
  km_by_grade,
  data         = GBSG2,
  pval         = TRUE,
  conf.int     = FALSE,
  legend.title = "Tumor grade",
  xlab         = "Time (days)",
  ylab         = "Recurrence-free survival probability",
  title        = "Kaplan-Meier Curves by Tumor Grade",
  ggtheme      = theme_minimal()
)
ggsave(
  filename = "outputs/figures/02_kaplan_meier_by_grade.png",
  plot     = grade_plot$plot,
  width    = 7, height = 6, dpi = 300
)
cat("Saved outputs/figures/02_kaplan_meier_by_grade.png\n")

# ---- 7. Save numeric results for the README / report -----------------------
km_results <- data.frame(
  comparison    = "Hormonal therapy vs. no therapy",
  chisq         = round(logrank_test$chisq, 3),
  df            = length(logrank_test$n) - 1,
  p_value       = round(logrank_p, 4)
)
write.csv(km_results, "outputs/tables/01_logrank_test_result.csv", row.names = FALSE)
print(km_results)
