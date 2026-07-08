############################################################
# analyse_hypotheses.R
# ----------------------------------------------------------
# Interprets pc_results in terms of the three simulation
# hypotheses (H1, H2, H3).
#
# Assumes pc_results is already in the workspace.
# If not, load it first:
#   pc_results = readRDS("pc_results_2500sim_....rds")
############################################################

options(scipen = 999)
n = nrow(pc_results)

cat("========================================================\n")
cat("SIMULATION RESULTS — HYPOTHESIS EVALUATION\n")
cat(sprintf("  n = %d iterations\n", n))
cat("========================================================\n\n")


# ==============================================================
# H1: Project-level UB coverage
# ==============================================================
# UB_coverage is the fraction of projects per iteration where
# theta_j <= UB_j. The nominal level is 0.95.
# H1,0: E[UB_coverage] = 0.95
# H1,1: E[UB_coverage] < 0.95  (one-sided)

cat("--------------------------------------------------------\n")
cat("H1: Project-level UB coverage (nominal = 0.95)\n")
cat("--------------------------------------------------------\n")

mean_cov  = mean(pc_results$UB_coverage, na.rm = TRUE)
sd_cov    = sd(pc_results$UB_coverage,   na.rm = TRUE)
se_cov    = sd_cov / sqrt(n)

# One-sample one-sided t-test: H1,1 is coverage < 0.95
t_h1      = (mean_cov - 0.95) / se_cov
p_h1      = pt(t_h1, df = n - 1)   # lower tail

cat(sprintf("  Mean UB_coverage : %.4f  (SD = %.4f)\n", mean_cov, sd_cov))
cat(sprintf("  t = %.3f,  one-sided p = %.4f\n", t_h1, p_h1))
cat(sprintf("  Proportion of iterations with coverage < 0.95: %.3f\n",
            mean(pc_results$UB_coverage < 0.95, na.rm = TRUE)))

if (p_h1 < 0.05) {
  cat("  -> H1,0 REJECTED: coverage is significantly below 0.95\n\n")
} else {
  cat("  -> H1,0 NOT rejected at alpha = 0.05\n\n")
}


# ==============================================================
# H2: Point-estimate accuracy at payment-claim level
# ==============================================================
# Bias: E[estimator] - theta_PC
# H2,0: E[theta_PC_UVB] = theta_PC  AND  E[theta_PC_ADR] = theta_PC
# H2,1: at least one is biased
# Additional: Var(ADR) > Var(UVB)

cat("--------------------------------------------------------\n")
cat("H2: Point-estimate accuracy\n")
cat("--------------------------------------------------------\n")

bias_uvb  = pc_results$theta_PC_UVB - pc_results$TER_true
bias_adr  = pc_results$theta_PC_ADR - pc_results$TER_true

# Two-sided t-tests for bias = 0
t_uvb     = t.test(bias_uvb, mu = 0)
t_adr     = t.test(bias_adr, mu = 0)

cat("  UVB estimator:\n")
cat(sprintf("    Mean bias : %+.5f  (SD = %.5f)\n",
            mean(bias_uvb, na.rm = TRUE), sd(bias_uvb, na.rm = TRUE)))
cat(sprintf("    t = %.3f,  two-sided p = %.4f\n",
            t_uvb$statistic, t_uvb$p.value))

cat("  ADR estimator:\n")
cat(sprintf("    Mean bias : %+.5f  (SD = %.5f)\n",
            mean(bias_adr, na.rm = TRUE), sd(bias_adr, na.rm = TRUE)))
cat(sprintf("    t = %.3f,  two-sided p = %.4f\n",
            t_adr$statistic, t_adr$p.value))

# Variance comparison (one-sided F-test: Var_ADR > Var_UVB)
var_uvb   = var(bias_uvb, na.rm = TRUE)
var_adr   = var(bias_adr, na.rm = TRUE)
f_stat    = var_adr / var_uvb
p_ftest   = pf(f_stat, df1 = n - 1, df2 = n - 1, lower.tail = FALSE)

cat(sprintf("\n  Var(UVB bias) = %.6f\n", var_uvb))
cat(sprintf("  Var(ADR bias) = %.6f\n", var_adr))
cat(sprintf("  F = %.3f,  one-sided p (ADR > UVB) = %.4f\n", f_stat, p_ftest))

if (p_ftest < 0.05) {
  cat("  -> Additional hypothesis SUPPORTED: Var(ADR) > Var(UVB)\n\n")
} else {
  cat("  -> Var(ADR) > Var(UVB) NOT supported at alpha = 0.05\n\n")
}


# ==============================================================
# H3: Operational consequences — unjustified adverse opinions
# ==============================================================
# Residual error not corrected by UVB: theta_PC - theta_PC_UVB
# Unjustified adverse opinion (RTER-based):
#   RTER > 0.02  AND  (theta_PC - theta_PC_UVB) <= 0.02
# Same definition with TER for comparison.
#
# H3,0: P(unjustified adverse | RTER) > 0   <- this is the null (expected positive)
# H3,1: P(unjustified adverse | RTER) <= 0  <- reject if never occurs

cat("--------------------------------------------------------\n")
cat("H3: Operational consequences — unjustified adverse opinions\n")
cat("--------------------------------------------------------\n")

residual = pc_results$TER_true - pc_results$theta_PC_UVB

adverse_rter   = pc_results$RTER > 0.02
adverse_ter    = pc_results$TER  > 0.02
unjustified    = residual <= 0.02

rate_rter_unjust = mean(adverse_rter & unjustified, na.rm = TRUE)
rate_ter_unjust  = mean(adverse_ter  & unjustified, na.rm = TRUE)
rate_rter_any    = mean(adverse_rter, na.rm = TRUE)
rate_ter_any     = mean(adverse_ter,  na.rm = TRUE)

cat(sprintf("  RTER > 2%%: %.3f of iterations\n", rate_rter_any))
cat(sprintf("    of which unjustified (residual <= 2%%): %.3f\n", rate_rter_unjust))
cat(sprintf("    absolute count: %d / %d\n",
            sum(adverse_rter & unjustified, na.rm = TRUE), n))

cat(sprintf("\n  TER > 2%%:  %.3f of iterations\n", rate_ter_any))
cat(sprintf("    of which unjustified (residual <= 2%%): %.3f\n", rate_ter_unjust))
cat(sprintf("    absolute count: %d / %d\n",
            sum(adverse_ter & unjustified, na.rm = TRUE), n))

cat(sprintf("\n  RTER unjustified rate <= TER unjustified rate: %s\n",
            if (rate_rter_unjust <= rate_ter_unjust) "YES" else "NO"))

if (rate_rter_unjust > 0) {
  cat("  -> H3,0 SUPPORTED: unjustified adverse rate (RTER) > 0\n\n")
} else {
  cat("  -> H3,0 NOT supported: no unjustified adverse opinions observed\n\n")
}

cat("========================================================\n")
cat("END OF HYPOTHESIS EVALUATION\n")
cat("========================================================\n")


# ==============================================================
# SECTION D — DESCRIPTIVE STATISTICS AND PLOTS
# ==============================================================
# All quantities are expressed as percentage points (x100)
# so axes read as "% of book value" rather than decimals.

to_pct = function(x) x * 100

bias_uvb_pct = to_pct(bias_uvb)   # theta_PC_UVB - TER_true
bias_adr_pct = to_pct(bias_adr)   # theta_PC_ADR - TER_true
residual_pct = to_pct(residual)   # TER_true - theta_PC_UVB  (UVB shortfall)

# Helper: compact quantile summary
quant_summary = function(x, label) {
  q = quantile(x, c(0.01, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99), na.rm = TRUE)
  cat(sprintf(
    "  %s\n    p01=%.2f  p05=%.2f  p25=%.2f  median=%.2f  p75=%.2f  p95=%.2f  p99=%.2f\n",
    label, q[1], q[2], q[3], q[4], q[5], q[6], q[7]
  ))
}

cat("\n--------------------------------------------------------\n")
cat("D1: UVB deviation from true TER  (theta_PC_UVB - TER_true, pct pts)\n")
cat("--------------------------------------------------------\n")
cat(sprintf("  Mean  : %+.3f pp\n", mean(bias_uvb_pct, na.rm = TRUE)))
cat(sprintf("  SD    : %.3f pp\n",  sd(bias_uvb_pct,   na.rm = TRUE)))
cat(sprintf("  RMSE  : %.3f pp\n",  sqrt(mean(bias_uvb_pct^2, na.rm = TRUE))))
cat(sprintf("  UVB underestimates (bias < 0): %.1f%%\n",
            100 * mean(bias_uvb_pct < 0, na.rm = TRUE)))
cat(sprintf("  UVB overestimates  (bias > 0): %.1f%%\n",
            100 * mean(bias_uvb_pct > 0, na.rm = TRUE)))
quant_summary(bias_uvb_pct, "Quantiles:")

cat("\n--------------------------------------------------------\n")
cat("D2: ADR deviation from true TER  (theta_PC_ADR - TER_true, pct pts)\n")
cat("--------------------------------------------------------\n")
cat(sprintf("  Mean  : %+.3f pp\n", mean(bias_adr_pct, na.rm = TRUE)))
cat(sprintf("  SD    : %.3f pp\n",  sd(bias_adr_pct,   na.rm = TRUE)))
cat(sprintf("  RMSE  : %.3f pp\n",  sqrt(mean(bias_adr_pct^2, na.rm = TRUE))))
cat(sprintf("  ADR above true error  (bias > 0): %.1f%%  of iterations\n",
            100 * mean(bias_adr_pct > 0, na.rm = TRUE)))
cat(sprintf("  ADR >2pp above true   (bias > 2): %.1f%%  of iterations\n",
            100 * mean(bias_adr_pct > 2, na.rm = TRUE)))
cat(sprintf("  ADR >5pp above true   (bias > 5): %.1f%%  of iterations\n",
            100 * mean(bias_adr_pct > 5, na.rm = TRUE)))

# Conditional magnitude: how far above when ADR > true?
adr_overest = bias_adr_pct[bias_adr_pct > 0]
cat(sprintf("  When ADR > true: mean excess = %.3f pp  (n = %d)\n",
            mean(adr_overest, na.rm = TRUE), length(adr_overest)))
quant_summary(bias_adr_pct, "Quantiles:")

cat("\n--------------------------------------------------------\n")
cat("D3: UVB residual shortfall  (TER_true - theta_PC_UVB, pct pts)\n")
cat("    Positive = UVB undercorrected; negative = UVB overcorrected\n")
cat("--------------------------------------------------------\n")
cat(sprintf("  Mean shortfall  : %+.3f pp\n", mean(residual_pct, na.rm = TRUE)))
cat(sprintf("  Shortfall > 0   : %.1f%%  of iterations\n",
            100 * mean(residual_pct > 0, na.rm = TRUE)))
cat(sprintf("  Shortfall > 2pp : %.1f%%  of iterations\n",
            100 * mean(residual_pct > 2, na.rm = TRUE)))
quant_summary(residual_pct, "Quantiles:")


# ==============================================================
# PLOTS — 2x3 grid saved to PDF
# ==============================================================

pdf("hypothesis_plots.pdf", width = 12, height = 8)
par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))

# -- Plot 1: UB_coverage distribution -------------------------
hist(to_pct(pc_results$UB_coverage),
     breaks = 40,
     main   = "H1: UB coverage per iteration",
     xlab   = "% projects with theta_j <= UB_j",
     col    = "steelblue", border = "white")
abline(v = 95, col = "red", lwd = 2, lty = 2)
legend("topleft", legend = "Nominal 95%", col = "red", lty = 2, lwd = 2, bty = "n")

# -- Plot 2: UVB bias distribution ----------------------------
hist(bias_uvb_pct,
     breaks = 40,
     main   = "H2: UVB bias  (theta_hat_UVB - TER_true)",
     xlab   = "Bias (percentage points)",
     col    = "steelblue", border = "white")
abline(v = 0, col = "red", lwd = 2, lty = 2)
abline(v = mean(bias_uvb_pct, na.rm = TRUE), col = "orange", lwd = 2)
legend("topright",
       legend = c("Zero", "Mean"),
       col    = c("red", "orange"), lty = 2, lwd = 2, bty = "n")

# -- Plot 3: ADR bias distribution ----------------------------
hist(bias_adr_pct,
     breaks = 40,
     main   = "H2: ADR bias  (theta_hat_ADR - TER_true)",
     xlab   = "Bias (percentage points)",
     col    = "steelblue", border = "white")
abline(v = 0, col = "red", lwd = 2, lty = 2)
abline(v = mean(bias_adr_pct, na.rm = TRUE), col = "orange", lwd = 2)
legend("topright",
       legend = c("Zero", "Mean"),
       col    = c("red", "orange"), lty = 2, lwd = 2, bty = "n")

# -- Plot 4: UVB vs ADR bias side-by-side boxplot -------------
boxplot(bias_uvb_pct, bias_adr_pct,
        names  = c("UVB", "ADR"),
        main   = "H2: Bias distributions compared",
        ylab   = "Bias (percentage points)",
        col    = c("steelblue", "coral"),
        outline = TRUE)
abline(h = 0, col = "red", lwd = 2, lty = 2)

# -- Plot 5: UVB residual shortfall ---------------------------
hist(residual_pct,
     breaks = 40,
     main   = "H3: UVB residual shortfall\n(TER_true - theta_hat_UVB)",
     xlab   = "Shortfall (percentage points)",
     col    = "steelblue", border = "white")
abline(v = 0,  col = "red",    lwd = 2, lty = 2)
abline(v = 2,  col = "orange", lwd = 2, lty = 2)
legend("topright",
       legend = c("Zero", "Materiality (2%)"),
       col    = c("red", "orange"), lty = 2, lwd = 2, bty = "n")

# -- Plot 6: RTER vs residual shortfall scatter ---------------
# Colour: red = unjustified adverse (RTER>2% but shortfall<=2%)
col_vec = ifelse(adverse_rter & unjustified, "red",
                 ifelse(adverse_rter,               "orange",
                        "steelblue"))
plot(residual_pct, to_pct(pc_results$RTER),
     pch  = 16, cex = 0.4,
     col  = adjustcolor(col_vec, alpha.f = 0.5),
     main = "H3: RTER vs UVB shortfall",
     xlab = "UVB shortfall (pp)",
     ylab = "RTER (pp)")
abline(h = 2, col = "red",    lwd = 1.5, lty = 2)
abline(v = 2, col = "orange", lwd = 1.5, lty = 2)
legend("topright",
       legend = c("RTER<=2%", "RTER>2% justified", "RTER>2% unjustified"),
       col    = c("steelblue", "orange", "red"),
       pch    = 16, bty = "n", cex = 0.8)

dev.off()
cat("\nPlots saved to hypothesis_plots.pdf\n")