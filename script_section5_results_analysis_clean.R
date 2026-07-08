############################################################
# script_section5_results_analysis.R
# ----------------------------------------------------------
# Analysis script for Chapter 5 of the thesis.
#
# Reads pc_results from run_simulation_parallel.R and produces
# the values, tables, and figures referenced in §5.1 to §5.5.
#
# Structure:
#   Section 0 — Setup and helper functions
#   Section 1 — §5.2 Aim 1: Coverage of UVB's upper bound
#   Section 2 — §5.3 Aim 2: Point-estimate accuracy
#   Section 3 — §5.4 Aim 3: Operational consequences
#   Section 4 — §5.5 Decomposition (L0 to L5)
#
# Style:
#   = for assignment (not <-)
#   |> for piping (not %>%)
############################################################

library(ggplot2)
library(dplyr)
library(tidyr)
library(flextable)
library(officer)

options(scipen = 999)

# ============================================================
# SECTION 0 — Setup and helper functions
# ============================================================

# --- 0.1 Output directory -----------------------------------
cfg_output_dir = "output_chapter5"
if (!dir.exists(cfg_output_dir)) dir.create(cfg_output_dir)

# --- 0.2 Load pc_results ------------------------------------
# pc_results_*.rds is created by run_simulation_parallel.R and
# contains one row per simulation iteration.
df = pc_results
n_sim = nrow(df)

cat("============================================================\n")
cat("Loaded:", n_sim, "iterations,", ncol(df), "columns\n")
cat("Failed (any NA):", sum(!complete.cases(df)), "\n")
cat("============================================================\n\n")

# --- 0.3 Helper functions for MCSE --------------------------
# MCSE for a mean of i.i.d. observations:    sd(x) / sqrt(n)
# MCSE for an empirical SD:                  sd(x) / sqrt(2(n-1))
# MCSE for a proportion:                     sqrt(p(1-p)/n)

mcse_mean = function(x) sd(x) / sqrt(length(x))
mcse_sd   = function(x) sd(x) / sqrt(2 * (length(x) - 1))
mcse_prop = function(p, n) sqrt(p * (1 - p) / n)

# --- 0.4 Format helpers -------------------------------------
fmt_num = function(x, digits = 4) sprintf(paste0("%.", digits, "f"), x)
fmt_pct = function(x, digits = 2) sprintf(paste0("%.", digits, "f%%"), 100 * x)

# Format mean with MCSE in parentheses, e.g. "0.8883 (0.0003)"
fmt_with_mcse = function(value, mcse, digits = 4) {
  sprintf("%s (%s)", fmt_num(value, digits), fmt_num(mcse, digits))
}

# Standard flextable styling for thesis tables (Calibri 11, APA-like)
style_flextable = function(ft, caption_text) {
  ft |>
    set_caption(
      caption = as_paragraph(
        as_chunk("Table [xx]: ",
                 props = fp_text(italic = TRUE, font.family = "Calibri", font.size = 11)),
        as_chunk(caption_text,
                 props = fp_text(italic = TRUE, font.family = "Calibri", font.size = 11))
      ),
      fp_p = fp_par(text.align = "center")
    ) |>
    fontsize(size = 11, part = "all") |>
    font(fontname = "Calibri", part = "all") |>
    autofit() |>
    border_remove() |>
    hline_top(part = "header",  border = fp_border(width = 1.5)) |>
    hline_bottom(part = "header", border = fp_border(width = 0.75)) |>
    hline_bottom(part = "body",   border = fp_border(width = 1.5))
}

# Standard ggplot theme for thesis figures
theme_thesis = function() {
  theme_classic(base_family = "sans", base_size = 11) +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey85", linewidth = 0.3),
      axis.line          = element_line(linewidth = 0.4),
      axis.ticks         = element_line(linewidth = 0.4),
      strip.background   = element_blank(),
      strip.text         = element_text(face = "bold"),
      legend.position    = "bottom",
      plot.margin        = margin(5, 10, 5, 5)
    )
}

# ============================================================
# SECTION 1 — §5.2 Aim 1: Coverage of UVB's upper bound
# ============================================================
# Test:
#   H_1,0: P(theta_j <= UB_j) = 0.95
#   H_1,1: P(theta_j <= UB_j) < 0.95   (one-sided)
#
# Output:
#   - Mean coverage, MCSE, one-sided 95% upper bound, z-stat, p-value
#   - Decision against H_1,0
#   - Figure: histogram of per-iteration coverage rates

cat("============================================================\n")
cat("§5.2 Aim 1: Coverage of UVB's upper bound\n")
cat("============================================================\n")

# --- 1.1 Point estimate and MCSE ----------------------------
nominal  = 0.95
mean_cov = mean(df$UB_coverage)
mcse_cov = mcse_mean(df$UB_coverage)

# --- 1.2 One-sided test against H_1,0 -----------------------
# z = (mean - nominal) / MCSE
# Under H_1,0, z is approximately standard normal (n_sim large).
# H_1,1 is mu < 0.95, so the p-value is the lower tail.
z_cov = (mean_cov - nominal) / mcse_cov
p_cov = pnorm(z_cov)

# 95% one-sided upper bound (matched to the one-sided test)
upper_bound_cov = mean_cov + qnorm(0.95) * mcse_cov

# Decision
deviation_pp_cov = (mean_cov - nominal) * 100   # in percentage points
if (p_cov < 0.05) {
  decision_h1 = "REJECTED at 5% significance level"
} else {
  decision_h1 = "NOT REJECTED at 5% significance level"
}

# --- 1.3 Console output -------------------------------------
cat("Mean coverage:              ", fmt_num(mean_cov, 4), "\n")
cat("MCSE coverage:              ", fmt_num(mcse_cov, 4), "\n")
cat("Nominal level:              ", nominal, "\n")
cat("Deviation (pp):             ", sprintf("%+.2f", deviation_pp_cov), "\n")
cat("One-sided 95% upper bound:  ", fmt_num(upper_bound_cov, 4), "\n")
cat("z-statistic:                ", fmt_num(z_cov, 1), "\n")
cat("p-value:                    ", ifelse(p_cov < 0.001, "< 0.001", fmt_num(p_cov, 4)), "\n")
cat("H_1,0 decision:             ", decision_h1, "\n\n")

# --- 1.4 Figure: histogram of per-iteration coverage --------
p_5_2 = ggplot(df, aes(x = UB_coverage)) +
  geom_histogram(binwidth = 0.01, boundary = 0,
                 fill = "#4393C3", color = "white", linewidth = 0.3) +
  geom_vline(xintercept = nominal,
             color = "darkred", linewidth = 0.8, linetype = "dashed") +
  geom_vline(xintercept = mean_cov,
             color = "black", linewidth = 0.5) +
  annotate("text", x = nominal,  y = Inf,
           label = "Nominal 0.95",
           hjust = -0.1, vjust = 1.5,
           color = "darkred", family = "sans", size = 3.5) +
  annotate("text", x = mean_cov, y = Inf,
           label = paste0("Mean ", fmt_num(mean_cov, 3)),
           hjust = 1.1,  vjust = 1.5,
           color = "black",  family = "sans", size = 3.5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = expression("Per-iteration coverage rate " * Coverage[i]),
       y = "Count") +
  theme_thesis()

ggsave(file.path(cfg_output_dir, "figure_5_2_coverage.png"),
       p_5_2, width = 8, height = 5, dpi = 600)
print(p_5_2)

# --- 1.5 Store results --------------------------------------
results_5_2 = list(
  mean_cov         = mean_cov,
  mcse_cov         = mcse_cov,
  upper_bound      = upper_bound_cov,
  deviation_pp     = deviation_pp_cov,
  z_stat           = z_cov,
  p_value          = p_cov,
  decision         = decision_h1
)

# ============================================================
# SECTION 2 — §5.3 Aim 2: Point-estimate accuracy
# ============================================================
# Tests:
#   H_2,0: E[theta_PC^UVB] = theta_PC AND E[theta_PC^ADR] = theta_PC
#   H_2,1: at least one estimator is biased  (two-sided per estimator)
# Additional descriptive comparison:
#   Var(theta_PC^ADR) > Var(theta_PC^UVB)  (no formal test)
#
# Output:
#   - Bias, EmpSE, RMSE for UVB and ADR with MCSE
#   - 95% two-sided MC CI for the bias of each estimator
#   - p-values and decision against H_2,0
#   - Variance comparison
#   - Table for §5.3
#   - Figure: density of deviations from truth

cat("============================================================\n")
cat("§5.3 Aim 2: Point-estimate accuracy\n")
cat("============================================================\n")

# --- 2.1 Per-iteration deviations from truth ---------------
df$dev_uvb = df$theta_PC_UVB - df$TER_true   # negative -> UVB underestimates
df$dev_adr = df$theta_PC_ADR - df$TER_true

# --- 2.2 Bias, EmpSE, RMSE ---------------------------------
# Bias:  mean of deviations
# EmpSE: SD of estimator across iterations (around its own mean)
# RMSE:  sqrt(mean of squared deviations)
bias_uvb = mean(df$dev_uvb)
bias_adr = mean(df$dev_adr)

empse_uvb = sd(df$theta_PC_UVB)
empse_adr = sd(df$theta_PC_ADR)

rmse_uvb = sqrt(mean(df$dev_uvb^2))
rmse_adr = sqrt(mean(df$dev_adr^2))

# --- 2.3 MCSEs ---------------------------------------------
mcse_bias_uvb  = empse_uvb / sqrt(n_sim)
mcse_bias_adr  = empse_adr / sqrt(n_sim)
mcse_empse_uvb = empse_uvb / sqrt(2 * (n_sim - 1))
mcse_empse_adr = empse_adr / sqrt(2 * (n_sim - 1))

# --- 2.4 Two-sided test of bias against zero ---------------
# z = bias / MCSE(bias);  p = 2 * pnorm(-|z|)
z_bias_uvb = bias_uvb / mcse_bias_uvb
z_bias_adr = bias_adr / mcse_bias_adr
p_bias_uvb = 2 * pnorm(-abs(z_bias_uvb))
p_bias_adr = 2 * pnorm(-abs(z_bias_adr))

# Two-sided 95% MC CI for the bias
ci_bias_uvb = c(bias_uvb - qnorm(0.975) * mcse_bias_uvb,
                bias_uvb + qnorm(0.975) * mcse_bias_uvb)
ci_bias_adr = c(bias_adr - qnorm(0.975) * mcse_bias_adr,
                bias_adr + qnorm(0.975) * mcse_bias_adr)

# H_2,0 is rejected if EITHER estimator has p < 0.05
reject_uvb = p_bias_uvb < 0.05
reject_adr = p_bias_adr < 0.05
if (reject_uvb || reject_adr) {
  decision_h2 = "REJECTED at 5% significance level (at least one estimator biased)"
} else {
  decision_h2 = "NOT REJECTED at 5% significance level"
}

# --- 2.5 Variance comparison (descriptive, no formal test) -
var_uvb   = empse_uvb^2
var_adr   = empse_adr^2
var_ratio = var_adr / var_uvb
if (var_adr > var_uvb) {
  var_verdict = "CONFIRMED (Var(ADR) > Var(UVB))"
} else {
  var_verdict = "REFUTED (Var(ADR) <= Var(UVB))"
}

# --- 2.6 Console output ------------------------------------
cat("\nUVB estimator:\n")
cat("  Bias:                ", fmt_num(bias_uvb, 5),
    "  (MCSE", fmt_num(mcse_bias_uvb, 5), ")\n")
cat("  EmpSE:               ", fmt_num(empse_uvb, 5),
    "  (MCSE", fmt_num(mcse_empse_uvb, 5), ")\n")
cat("  RMSE:                ", fmt_num(rmse_uvb, 5), "\n")
cat("  95% MC CI for bias:  [", fmt_num(ci_bias_uvb[1], 5), ",",
    fmt_num(ci_bias_uvb[2], 5), "]\n")
cat("  z-statistic:         ", fmt_num(z_bias_uvb, 3), "\n")
cat("  p-value (two-sided): ", ifelse(p_bias_uvb < 0.001, "< 0.001",
                                      fmt_num(p_bias_uvb, 4)), "\n")

cat("\nADR estimator:\n")
cat("  Bias:                ", fmt_num(bias_adr, 5),
    "  (MCSE", fmt_num(mcse_bias_adr, 5), ")\n")
cat("  EmpSE:               ", fmt_num(empse_adr, 5),
    "  (MCSE", fmt_num(mcse_empse_adr, 5), ")\n")
cat("  RMSE:                ", fmt_num(rmse_adr, 5), "\n")
cat("  95% MC CI for bias:  [", fmt_num(ci_bias_adr[1], 5), ",",
    fmt_num(ci_bias_adr[2], 5), "]\n")
cat("  z-statistic:         ", fmt_num(z_bias_adr, 3), "\n")
cat("  p-value (two-sided): ", ifelse(p_bias_adr < 0.001, "< 0.001",
                                      fmt_num(p_bias_adr, 4)), "\n")

cat("\nH_2,0 decision:        ", decision_h2, "\n")
cat("\nVariance comparison:\n")
cat("  Var(UVB):            ", fmt_num(var_uvb, 6), "\n")
cat("  Var(ADR):            ", fmt_num(var_adr, 6), "\n")
cat("  Ratio Var(ADR)/Var(UVB):", fmt_num(var_ratio, 3), "\n")
cat("  Verdict:             ", var_verdict, "\n\n")

# --- 2.7 Table for §5.3 ------------------------------------
table_5_3 = data.frame(
  Estimator = c("UVB", "ADR"),
  Bias      = c(fmt_with_mcse(bias_uvb,  mcse_bias_uvb),
                fmt_with_mcse(bias_adr,  mcse_bias_adr)),
  EmpSE     = c(fmt_with_mcse(empse_uvb, mcse_empse_uvb),
                fmt_with_mcse(empse_adr, mcse_empse_adr)),
  RMSE      = c(fmt_num(rmse_uvb, 4), fmt_num(rmse_adr, 4)),
  check.names = FALSE
)
names(table_5_3) = c("Estimator", "Bias (MCSE)", "EmpSE (MCSE)", "RMSE")

ft_5_3 = flextable(table_5_3) |>
  align(j = 1,    align = "left",  part = "all") |>
  align(j = 2:4,  align = "right", part = "all") |>
  style_flextable(
    "Bias, empirical standard error, and root mean squared error of UVB and ADR point estimates relative to the true payment-claim error rate. Monte Carlo standard errors in parentheses."
  )

save_as_docx(ft_5_3, path = file.path(cfg_output_dir, "table_5_3_accuracy.docx"))
ft_5_3

# --- 2.8 Figure: density of deviations from truth ----------
# Plotting deviations rather than absolute estimates centres
# both methods on zero, isolating sampling error from the
# variation in theta_PC across iterations.
plot_dev = data.frame(
  deviation = c(df$dev_uvb, df$dev_adr),
  method    = factor(rep(c("UVB", "ADR"), each = n_sim),
                     levels = c("UVB", "ADR"))
)

p_5_3 = ggplot(plot_dev, aes(x = deviation, fill = method)) +
  geom_density(alpha = 0.5, color = NA) +
  geom_vline(xintercept = 0, color = "black",
             linewidth = 0.5, linetype = "dashed") +
  scale_fill_manual(values = c("UVB" = "#4393C3", "ADR" = "#C0392B")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = expression(hat(theta)[PC] - theta[PC]),
       y = "Density",
       fill = NULL) +
  theme_thesis()

ggsave(file.path(cfg_output_dir, "figure_5_3_accuracy.png"),
       p_5_3, width = 8, height = 5, dpi = 600)
print(p_5_3)

# --- 2.9 Store results -------------------------------------
results_5_3 = list(
  bias_uvb  = bias_uvb,  mcse_bias_uvb  = mcse_bias_uvb,
  bias_adr  = bias_adr,  mcse_bias_adr  = mcse_bias_adr,
  empse_uvb = empse_uvb, mcse_empse_uvb = mcse_empse_uvb,
  empse_adr = empse_adr, mcse_empse_adr = mcse_empse_adr,
  rmse_uvb  = rmse_uvb,  rmse_adr       = rmse_adr,
  ci_bias_uvb = ci_bias_uvb,
  ci_bias_adr = ci_bias_adr,
  z_bias_uvb  = z_bias_uvb,  p_bias_uvb = p_bias_uvb,
  z_bias_adr  = z_bias_adr,  p_bias_adr = p_bias_adr,
  decision    = decision_h2,
  var_ratio   = var_ratio,   var_verdict = var_verdict
)

# ============================================================
# SECTION 3 — §5.4 Aim 3: Operational consequences
# ============================================================
# Test:
#   H_3,0: P(RTER > 0.02 AND theta_PC - theta_PC^UVB <= 0.02) = 0
#   H_3,1: P(RTER > 0.02 AND theta_PC - theta_PC^UVB <= 0.02) > 0
#
# Notes on definitions (matches §4.3 of the thesis):
#   - Justified adverse:   RTER > 0.02  AND  theta_PC - theta_UVB > 0.02
#                          (UVB undershot truth by more than 2 pp)
#   - Unjustified adverse: RTER > 0.02  AND  theta_PC - theta_UVB <= 0.02
#                          (UVB did NOT undershoot by more than 2 pp)
#
# In the data, df$diff has the opposite sign convention:
#   diff = theta_PC_UVB - TER_true
# so undershoot > 2pp corresponds to diff < -0.02.

cat("============================================================\n")
cat("§5.4 Aim 3: Operational consequences\n")
cat("============================================================\n")

cfg_threshold = 0.02

# --- 3.1 Indicator construction ----------------------------
df$diff_uvb_truth = df$theta_PC_UVB - df$TER_true   # < -0.02 means UVB undershot

adverse_rter   = df$RTER > cfg_threshold
uvb_undershoot = df$diff_uvb_truth < -cfg_threshold

just_rter   = adverse_rter & uvb_undershoot
unjust_rter = adverse_rter & !uvb_undershoot

# --- 3.2 Proportions ---------------------------------------
ar_rter        = mean(adverse_rter)
just_ar_rter   = mean(just_rter)
unjust_ar_rter = mean(unjust_rter)

# --- 3.3 MCSEs for proportions -----------------------------
mcse_ar_rter        = mcse_prop(ar_rter,        n_sim)
mcse_just_ar_rter   = mcse_prop(just_ar_rter,   n_sim)
mcse_unjust_ar_rter = mcse_prop(unjust_ar_rter, n_sim)

# --- 3.4 One-sided test for H_3,0 --------------------------
# z = (p_hat - 0) / MCSE,  p-value = upper tail (H_3,1: p > 0)
z_h3 = unjust_ar_rter / mcse_unjust_ar_rter
p_h3 = 1 - pnorm(z_h3)

# Reject H_3,0 if 0 lies below the 95% one-sided lower bound
unjust_lower_bound = max(0, unjust_ar_rter - qnorm(0.95) * mcse_unjust_ar_rter)
n_unjust_rter      = sum(unjust_rter)

if (p_h3 < 0.05) {
  decision_h3 = "REJECTED at 5% significance level"
} else {
  decision_h3 = "NOT REJECTED at 5% significance level"
}

# --- 3.5 Console output ------------------------------------
cat("\nRTER (RTER > 2%):\n")
cat("  Adverse rate:               ", fmt_pct(ar_rter),
    "  (MCSE", fmt_num(mcse_ar_rter, 4), ")\n")
cat("  Justified adverse rate:     ", fmt_pct(just_ar_rter),
    "  (MCSE", fmt_num(mcse_just_ar_rter, 4), ")\n")
cat("  Unjustified adverse rate:   ", fmt_pct(unjust_ar_rter),
    "  (MCSE", fmt_num(mcse_unjust_ar_rter, 4), ")\n")
cat("  Number of unjustified cases:", n_unjust_rter, "/", n_sim, "\n")
cat("  One-sided 95% lower bound:  ", fmt_pct(unjust_lower_bound), "\n")
cat("  H_3,0 decision:             ", decision_h3, "\n\n")
cat("  z-statistic:                ", fmt_num(z_h3, 1), "\n")
cat("  p-value:                    ", ifelse(p_h3 < 0.001, "< 0.001", fmt_num(p_h3, 4)), "\n")

# --- 3.6 Table for §5.4 ------------------------------------
table_5_4 = data.frame(
  Quantity = c("Adverse rate",
               "Justified adverse rate",
               "Unjustified adverse rate"),
  Value    = c(fmt_with_mcse(ar_rter,        mcse_ar_rter,        digits = 4),
               fmt_with_mcse(just_ar_rter,   mcse_just_ar_rter,   digits = 4),
               fmt_with_mcse(unjust_ar_rter, mcse_unjust_ar_rter, digits = 4)),
  check.names = FALSE
)
names(table_5_4) = c("Quantity", "Value (MCSE)")

ft_5_4 = flextable(table_5_4) |>
  align(j = 1, align = "left",  part = "all") |>
  align(j = 2, align = "right", part = "all") |>
  style_flextable(
    "Adverse rate, justified adverse rate, and unjustified adverse rate against the RTER. Monte Carlo standard errors in parentheses."
  )

save_as_docx(ft_5_4, path = file.path(cfg_output_dir, "table_5_4_operational.docx"))
ft_5_4

# --- 3.7 Store results -------------------------------------
results_5_4 = list(
  ar_rter            = ar_rter,        mcse_ar_rter        = mcse_ar_rter,
  just_ar_rter       = just_ar_rter,   mcse_just_ar_rter   = mcse_just_ar_rter,
  unjust_ar_rter     = unjust_ar_rter, mcse_unjust_ar_rter = mcse_unjust_ar_rter,
  n_unjust_rter      = n_unjust_rter,
  unjust_lower_bound = unjust_lower_bound,
  z_h3               = z_h3,
  p_h3               = p_h3,
  decision_h3        = decision_h3
)

# ============================================================
# SECTION 4 — §5.5 Decomposition (L0 to L5, residual scale)
# ============================================================
# Five levels on the residual scale (what UVB has not corrected):
#   L0_new = TER_true - theta_PC_UVB                   (truth)
#   L1_new = L1 - extrap_EE_UVB                        (stage-1 selection + extrapolation)
#   L2_new = L3                                         (+ within-project audit by ADR)
#   L3_new = L4                                         (+ incremental floor; = TER)
#   L4_new = L5                                         (+ denominator shrinkage; = RTER)
#
# All five levels are on the same residual scale, directly
# comparable to the 2% materiality threshold. Adverse rate at
# each level is therefore reported.
#
# Output:
#   - Console summary
#   - Table for §5.5
#   - Boxplot per level (with 2% threshold)
#   - Density plot of L0, L3, L4 (truth, TER, RTER)

cat("============================================================\n")
cat("§5.5 Decomposition (residual scale)\n")
cat("============================================================\n")

# --- 4.1 Construct the new residual-scale ladder -----------
df$L0_new = df$TER_true - df$theta_PC_UVB     # truth: residual UVB has not corrected
df$L1_new = df$L1 - df$extrap_EE_UVB          # stage-1 selection + extrapolation, on residual scale
df$L2_new = df$L3                             # = Extrap(EE_ADR - EE_UVB)
df$L3_new = df$L4                             # = TER
df$L4_new = df$L5                             # = RTER

levels_new = c("L0", "L1", "L2", "L3", "L4")  # display labels (no "_new" in output)
levels_col = c("L0_new", "L1_new", "L2_new", "L3_new", "L4_new")

# --- 4.2 Per-level summary statistics ----------------------
decomp = data.frame(
  Level       = levels_new,
  Mean        = sapply(levels_col, function(L) mean(df[[L]])),
  SD          = sapply(levels_col, function(L) sd(df[[L]])),
  MCSE        = sapply(levels_col, function(L) mcse_mean(df[[L]])),
  AdverseRate = sapply(levels_col, function(L) mean(df[[L]] > cfg_threshold)),
  row.names   = NULL
)
decomp$Delta_Mean = c(NA, diff(decomp$Mean))

# MCSE for adverse rate (a proportion)
decomp$MCSE_AR = sapply(decomp$AdverseRate, function(p) mcse_prop(p, n_sim))

# --- 4.3 Console output ------------------------------------
cat("\nDecomposition results (per level):\n")
for (i in seq_len(nrow(decomp))) {
  cat(sprintf(
    "  %s: mean = %s (MCSE %s), SD = %s, adverse rate = %s (MCSE %s)%s\n",
    decomp$Level[i],
    fmt_num(decomp$Mean[i], 5),
    fmt_num(decomp$MCSE[i], 5),
    fmt_num(decomp$SD[i], 5),
    fmt_pct(decomp$AdverseRate[i]),
    fmt_num(decomp$MCSE_AR[i], 4),
    if (i > 1) sprintf(" | delta-mean from previous = %+.5f", decomp$Delta_Mean[i]) else ""
  ))
}
cat("\n")

# --- 4.4 Table for §5.5 ------------------------------------
table_5_5 = data.frame(
  Level     = decomp$Level,
  Mean_MCSE = mapply(fmt_with_mcse, decomp$Mean, decomp$MCSE, digits = 4),
  SD        = fmt_num(decomp$SD, 4),
  Delta     = ifelse(is.na(decomp$Delta_Mean), "\u2014",
                     sprintf("%+s", fmt_num(decomp$Delta_Mean, 4))),
  Adverse   = mapply(function(p, m) sprintf("%s (%s)", fmt_pct(p), fmt_num(m, 4)),
                     decomp$AdverseRate, decomp$MCSE_AR),
  check.names = FALSE
)
names(table_5_5) = c("Level", "Mean (MCSE)", "SD",
                     "\u0394 Mean from previous", "Adverse rate (MCSE)")

ft_5_5 = flextable(table_5_5) |>
  align(j = 1,   align = "left",  part = "all") |>
  align(j = 2:5, align = "right", part = "all") |>
  style_flextable(
    "Decomposition of the path from the true UVB residual to the RTER. Mean, empirical standard deviation, Monte Carlo standard error, change in mean from the preceding level, and adverse rate (proportion of iterations exceeding the 2% materiality threshold) for each level L_0 to L_4."
  )

save_as_docx(ft_5_5, path = file.path(cfg_output_dir, "table_5_5_decomposition.docx"))
ft_5_5

# --- 4.5 Boxplot per level ---------------------------------
plot_decomp = data.frame(
  level = factor(rep(levels_new, each = n_sim), levels = levels_new),
  value = c(df$L0_new, df$L1_new, df$L2_new, df$L3_new, df$L4_new)
)

p_5_5 = ggplot(plot_decomp, aes(x = level, y = value)) +
  geom_boxplot(fill = "#4393C3", color = "black",
               outlier.size = 0.5, outlier.alpha = 0.3,
               linewidth = 0.4) +
  geom_hline(yintercept = cfg_threshold,
             color = "darkred", linewidth = 0.6, linetype = "dashed") +
  annotate("text", x = 0.7, y = cfg_threshold,
           label = "2% threshold",
           hjust = 0, vjust = -0.5,
           color = "darkred", family = "sans", size = 3.5) +
  labs(x = "Level", y = "Estimator value (residual scale)") +
  theme_thesis()

ggsave(file.path(cfg_output_dir, "figure_5_5_decomposition.png"),
       p_5_5, width = 9, height = 5, dpi = 600)
print(p_5_5)

# --- 4.6 Density plot: truth, TER, RTER --------------------
# Visualises the floor (transition L2 -> L3) and the
# denominator shrinkage (transition L3 -> L4) by overlaying
# the densities of L0 (truth), L3 (TER) and L4 (RTER).
plot_residual = data.frame(
  value    = c(df$L0_new, df$L3_new, df$L4_new),
  quantity = factor(rep(c("L0 (true UVB residual)",
                          "L3 (TER)",
                          "L4 (RTER)"),
                        each = n_sim),
                    levels = c("L0 (true UVB residual)",
                               "L3 (TER)",
                               "L4 (RTER)"))
)

p_5_5_density = ggplot(plot_residual,
                       aes(x = value, fill = quantity, color = quantity)) +
  geom_density(alpha = 0.3, linewidth = 0.6) +
  geom_vline(xintercept = cfg_threshold,
             color = "darkred", linewidth = 0.6, linetype = "dashed") +
  annotate("text", x = cfg_threshold, y = Inf,
           label = "2% threshold",
           hjust = -0.1, vjust = 1.5,
           color = "darkred", family = "sans", size = 3.5) +
  scale_fill_manual(values  = c("L0 (true UVB residual)" = "#2D6A4F",
                                "L3 (TER)"               = "#C0392B",
                                "L4 (RTER)"              = "#4393C3")) +
  scale_color_manual(values = c("L0 (true UVB residual)" = "#2D6A4F",
                                "L3 (TER)"               = "#C0392B",
                                "L4 (RTER)"              = "#4393C3")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Residual error rate",
       y = "Density",
       fill  = NULL,
       color = NULL) +
  theme_thesis()

ggsave(file.path(cfg_output_dir, "figure_5_5_density.png"),
       p_5_5_density, width = 9, height = 5, dpi = 600)
print(p_5_5_density)



# =============================================================================
# Figure 5.5b — Empirical density at each decomposition level (L0-L4)
# v4: per-panel blue label position; title/subtitle removed from figure
# =============================================================================

library(ggplot2)
library(patchwork)

windowsFonts(Calibri = windowsFont("Calibri"))

cfg_threshold = 0.02
cfg_output    = "figure_5_5_empirical_density_panels.png"

# --- Summary statistics ------------------------------------------------------

level_stats = data.frame(
  label = c("L0", "L1", "L2", "L3", "L4"),
  col   = c("L0_new", "L1_new", "L2_new", "L3_new", "L4_new")
)

level_stats$mu      = sapply(level_stats$col, function(v) mean(df[[v]]))
level_stats$sigma   = sapply(level_stats$col, function(v) sd(df[[v]]))
level_stats$p_above = sapply(level_stats$col, function(v) mean(df[[v]] >  cfg_threshold))
level_stats$p_below = sapply(level_stats$col, function(v) mean(df[[v]] <= cfg_threshold))

# -----------------------------------------------------------------------------
# BLUE LABEL POSITION — adjust x_lbl_left per panel here (as a fraction of
# the x-axis range, measured from the left edge to the threshold).
# 0.5 = halfway between left edge and threshold; higher = closer to threshold.
# -----------------------------------------------------------------------------
level_stats$x_lbl_left_frac = c(
  L0 = 0.45,   # adjust if label overlaps density peak
  L1 = 0.55,
  L2 = 0.68,
  L3 = 0.70,
  L4 = 0.57
)

# --- Plot function -----------------------------------------------------------

make_panel_empirical = function(label, values, p_above, p_below,
                                x_lbl_left_frac,
                                threshold = cfg_threshold) {
  
  kd     = density(values, n = 2048)
  df_kde = data.frame(x = kd$x, y = kd$y)
  
  y_at_threshold = approx(kd$x, kd$y, xout = threshold)$y
  
  df_left  = rbind(
    df_kde[df_kde$x <= threshold, ],
    data.frame(x = threshold, y = y_at_threshold)
  )
  df_right = rbind(
    data.frame(x = threshold, y = y_at_threshold),
    df_kde[df_kde$x >= threshold, ]
  )
  
  y_peak = max(kd$y)
  y_lbl  = 0.50 * y_peak
  sigma  = sd(values)
  x_lo   = min(kd$x)
  x_hi   = max(kd$x)
  
  # Blue label x: fraction of the distance from x_lo to threshold
  x_lbl_left  = x_lo + x_lbl_left_frac * (threshold - x_lo)
  
  # Red label x: threshold + 1.2 sigma, clamped
  x_lbl_right = min(threshold + 1.2 * sigma, x_hi - 0.05 * (x_hi - x_lo))
  
  ggplot() +
    geom_area(data = df_left,  aes(x = x, y = y),
              fill = "#4393C3", alpha = 0.25) +
    geom_area(data = df_right, aes(x = x, y = y),
              fill = "#D73027", alpha = 0.25) +
    geom_line(data = df_kde,   aes(x = x, y = y),
              linewidth = 0.6, color = "grey20") +
    geom_vline(xintercept = threshold,
               color = "red", linetype = "dashed", linewidth = 0.55) +
    annotate("text",
             x = x_lbl_left, y = y_lbl,
             label    = paste0(formatC(p_below * 100, format = "f", digits = 2), "%"),
             color    = "#2166AC", size = 3.2,
             family   = "Calibri", fontface = "bold") +
    annotate("text",
             x = x_lbl_right, y = y_lbl,
             label    = paste0(formatC(p_above * 100, format = "f", digits = 2), "%"),
             color    = "#B2182B", size = 3.2,
             family   = "Calibri", fontface = "bold") +
    scale_x_continuous(
      labels = scales::label_percent(accuracy = 0.1),
      expand = expansion(mult = 0.02)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(title = label, x = "Residual error rate", y = "Density") +
    theme_bw(base_family = "Calibri", base_size = 11) +
    theme(
      panel.border       = element_blank(),
      axis.line.x.bottom = element_line(linewidth = 0.4),
      axis.line.y.left   = element_line(linewidth = 0.4),
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.ticks         = element_line(linewidth = 0.4),
      strip.background   = element_blank(),
      plot.title         = element_text(face = "bold", size = 11),
      plot.margin        = margin(5, 10, 5, 5)
    )
}

# --- Build panels ------------------------------------------------------------

panels = vector("list", nrow(level_stats))

for (i in seq_len(nrow(level_stats))) {
  panels[[i]] = make_panel_empirical(
    label            = level_stats$label[i],
    values           = df[[level_stats$col[i]]],
    p_above          = level_stats$p_above[i],
    p_below          = level_stats$p_below[i],
    x_lbl_left_frac  = level_stats$x_lbl_left_frac[i]
  )
}

# --- Combine (no overall title/subtitle) -------------------------------------

combined = (panels[[1]] | panels[[2]] | panels[[3]]) /
  (panels[[4]] | panels[[5]] | plot_spacer())

# --- Save --------------------------------------------------------------------

ggsave(
  filename = cfg_output,
  plot     = combined,
  width    = 14,
  height   = 8,
  dpi      = 600
)

message("Saved: ", cfg_output)
combined

# --- 4.7 Store results -------------------------------------
results_5_5 = list(decomp = decomp)

# ============================================================
# END OF SCRIPT
# ============================================================
cat("============================================================\n")
cat("Analysis complete. All tables and figures saved to:\n")
cat("  ", normalizePath(cfg_output_dir), "\n")
cat("============================================================\n")
