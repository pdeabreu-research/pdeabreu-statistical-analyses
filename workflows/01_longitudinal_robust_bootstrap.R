# Longitudinal panel analysis with robust inference and bootstrap indirect effects
# Author: Pedro De Abreu
# Public portfolio version: sanitized and generalized from prior research workflows.
# Data: synthetic by default; no participant-level or proprietary data are included.

required <- c("sandwich", "lmtest", "boot")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop(
    "Missing required packages: ", paste(missing, collapse = ", "),
    "\nInstall with: install.packages(c(",
    paste(sprintf('"%s"', missing), collapse = ", "), "))"
  )
}

library(sandwich)
library(lmtest)
library(boot)

set.seed(20260823)
dir.create("outputs", showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# 1. Synthetic demonstration data
# -----------------------------------------------------------------------------
# Replace this block with a read.csv()/readr::read_csv() call for an authorized
# analytic dataset. Variable names are deliberately generic in this public demo.

simulate_panel_data <- function(n = 600L) {
  x_t1 <- rnorm(n)
  covariate_t1 <- rnorm(n)
  facet_a_t1 <- 0.65 * x_t1 + rnorm(n, sd = 0.75)
  facet_b_t1 <- 0.45 * x_t1 + rnorm(n, sd = 0.85)

  mediator_pos_t1 <- -0.30 * x_t1 + rnorm(n)
  mediator_neg_t1 <-  0.28 * x_t1 + rnorm(n)

  mediator_pos_t2 <- 0.58 * mediator_pos_t1 - 0.22 * x_t1 +
    0.10 * covariate_t1 + rnorm(n, sd = 0.80)
  mediator_neg_t2 <- 0.55 * mediator_neg_t1 + 0.20 * x_t1 +
    0.08 * covariate_t1 + rnorm(n, sd = 0.85)

  x_t2 <- 0.70 * x_t1 - 0.07 * mediator_pos_t1 +
    0.05 * mediator_neg_t1 + rnorm(n, sd = 0.70)

  outcome_help_t2 <- -0.12 * x_t1 + 0.30 * mediator_pos_t2 -
    0.14 * mediator_neg_t2 + 0.08 * covariate_t1 + rnorm(n)
  outcome_harm_t2 <- 0.14 * x_t1 - 0.24 * mediator_pos_t2 +
    0.28 * mediator_neg_t2 + 0.08 * covariate_t1 + rnorm(n)

  data.frame(
    participant_id = seq_len(n), x_t1, x_t2, covariate_t1,
    facet_a_t1, facet_b_t1,
    mediator_pos_t1, mediator_pos_t2,
    mediator_neg_t1, mediator_neg_t2,
    outcome_help_t2, outcome_harm_t2
  )
}

D <- simulate_panel_data()

# -----------------------------------------------------------------------------
# 2. Helpers: complete-case standardization and HC3 robust inference
# -----------------------------------------------------------------------------

z <- function(x) as.numeric(scale(x))

standardize_complete <- function(data, vars) {
  d <- data[complete.cases(data[vars]), vars, drop = FALSE]
  d[] <- lapply(d, z)
  d
}

fit_hc3 <- function(formula, data = D) {
  vars <- all.vars(formula)
  d <- standardize_complete(data, vars)
  fit <- lm(formula, data = d)
  V <- vcovHC(fit, type = "HC3")
  list(
    fit = fit,
    coefficients = coeftest(fit, vcov. = V),
    confint = coefci(fit, vcov. = V),
    n = nrow(d),
    r_squared = summary(fit)$r.squared
  )
}

# -----------------------------------------------------------------------------
# 3. Baseline-adjusted prospective and reciprocal-direction models
# -----------------------------------------------------------------------------

models <- list(
  positive_path = fit_hc3(mediator_pos_t2 ~ mediator_pos_t1 + x_t1 + covariate_t1),
  negative_path = fit_hc3(mediator_neg_t2 ~ mediator_neg_t1 + x_t1 + covariate_t1),
  reciprocal_x_from_positive = fit_hc3(x_t2 ~ x_t1 + mediator_pos_t1 + covariate_t1),
  reciprocal_x_from_negative = fit_hc3(x_t2 ~ x_t1 + mediator_neg_t1 + covariate_t1),
  help_outcome = fit_hc3(
    outcome_help_t2 ~ x_t1 + mediator_pos_t1 + mediator_pos_t2 +
      mediator_neg_t1 + mediator_neg_t2 + covariate_t1
  ),
  harm_outcome = fit_hc3(
    outcome_harm_t2 ~ x_t1 + mediator_pos_t1 + mediator_pos_t2 +
      mediator_neg_t1 + mediator_neg_t2 + covariate_t1
  ),
  facet_test = fit_hc3(mediator_pos_t2 ~ mediator_pos_t1 + facet_a_t1 + facet_b_t1)
)

model_table <- do.call(rbind, lapply(names(models), function(name) {
  x <- models[[name]]
  out <- as.data.frame(x$coefficients)
  out$term <- rownames(out)
  rownames(out) <- NULL
  out$model <- name
  out$n <- x$n
  out$r_squared <- x$r_squared
  out[, c("model", "term", "Estimate", "Std. Error", "t value", "Pr(>|t|)", "n", "r_squared")]
}))
write.csv(model_table, "outputs/longitudinal_hc3_models.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# 4. Parallel indirect effects with participant-level bootstrap
# -----------------------------------------------------------------------------
# These are prospective indirect associations: the mediator and outcome are
# both measured at Time 2. The workflow does not label them causal mediation.

parallel_indirect_stat <- function(data, indices, outcome) {
  d <- data[indices, , drop = FALSE]
  needed <- c(
    "x_t1", "covariate_t1",
    "mediator_pos_t1", "mediator_pos_t2",
    "mediator_neg_t1", "mediator_neg_t2",
    outcome
  )
  d <- standardize_complete(d, needed)

  m_pos <- lm(mediator_pos_t2 ~ mediator_pos_t1 + x_t1 + covariate_t1, data = d)
  m_neg <- lm(mediator_neg_t2 ~ mediator_neg_t1 + x_t1 + covariate_t1, data = d)
  y_fit <- lm(
    as.formula(paste(
      outcome,
      "~ x_t1 + mediator_pos_t1 + mediator_pos_t2 +",
      "mediator_neg_t1 + mediator_neg_t2 + covariate_t1"
    )),
    data = d
  )

  a_pos <- coef(m_pos)["x_t1"]
  a_neg <- coef(m_neg)["x_t1"]
  b_pos <- coef(y_fit)["mediator_pos_t2"]
  b_neg <- coef(y_fit)["mediator_neg_t2"]
  direct <- coef(y_fit)["x_t1"]

  c(
    a_pos = a_pos,
    b_pos = b_pos,
    indirect_pos = a_pos * b_pos,
    a_neg = a_neg,
    b_neg = b_neg,
    indirect_neg = a_neg * b_neg,
    indirect_difference = (a_pos * b_pos) - (a_neg * b_neg),
    direct = direct
  )
}

summarize_boot <- function(b) {
  estimates <- b$t0
  ci95 <- t(vapply(
    seq_along(estimates),
    function(j) quantile(b$t[, j], c(.025, .975), na.rm = TRUE),
    numeric(2)
  ))
  p_boot <- vapply(seq_along(estimates), function(j) {
    x <- b$t[, j]
    min(1, 2 * min(mean(x <= 0, na.rm = TRUE), mean(x >= 0, na.rm = TRUE)))
  }, numeric(1))

  data.frame(
    term = names(estimates),
    estimate = unname(estimates),
    ci95_low = ci95[, 1],
    ci95_high = ci95[, 2],
    bootstrap_p = p_boot,
    row.names = NULL
  )
}

R_BOOT <- 5000L
bootstrap_results <- list()
for (outcome in c("outcome_help_t2", "outcome_harm_t2")) {
  needed <- c(
    "x_t1", "covariate_t1",
    "mediator_pos_t1", "mediator_pos_t2",
    "mediator_neg_t1", "mediator_neg_t2",
    outcome
  )
  analytic <- D[complete.cases(D[needed]), , drop = FALSE]
  b <- boot(
    analytic,
    statistic = function(data, i) parallel_indirect_stat(data, i, outcome),
    R = R_BOOT
  )
  tmp <- summarize_boot(b)
  tmp$outcome <- outcome
  tmp$n <- nrow(analytic)
  bootstrap_results[[outcome]] <- tmp
}

write.csv(
  do.call(rbind, bootstrap_results),
  "outputs/longitudinal_bootstrap_indirect_effects.csv",
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 5. Exploratory correlation screen with Benjamini-Hochberg FDR control
# -----------------------------------------------------------------------------

predictors <- c("x_t1", "facet_a_t1", "facet_b_t1", "mediator_pos_t1", "mediator_neg_t1")
outcomes <- c("x_t2", "mediator_pos_t2", "mediator_neg_t2", "outcome_help_t2", "outcome_harm_t2")

screen <- do.call(rbind, lapply(predictors, function(pred) {
  do.call(rbind, lapply(outcomes, function(out) {
    ok <- complete.cases(D[[pred]], D[[out]])
    tst <- cor.test(D[[pred]][ok], D[[out]][ok])
    data.frame(
      predictor = pred,
      outcome = out,
      n = sum(ok),
      r = unname(tst$estimate),
      p = tst$p.value
    )
  }))
}))
screen$fdr_q <- p.adjust(screen$p, method = "BH")
screen$fdr_significant <- screen$fdr_q < .05
write.csv(screen, "outputs/exploratory_fdr_screen.csv", row.names = FALSE)

message("Complete. Outputs written to ./outputs/")
