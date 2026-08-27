# Reproducible example: time-dependent AUC analysis in two cohorts
#
# This example compares:
#   1. Biomarker alone
#   2. Clinical risk group alone
#   3. Biomarker + clinical risk group
#
# The AUC comparison tests clinical risk alone versus the combined model.

library(survival)
library(timeROC)
library(ggplot2)
library(dplyr)

set.seed(2026)

# ==============================================================================
# 1. Example data
# ==============================================================================

simulate_cohort <- function(n) {
  biomarker_score <- rnorm(n)
  clinical_risk_group <- sample(1:4, n, replace = TRUE)

  event_rate <- exp(
    -7.0 +
      0.35 * biomarker_score +
      0.25 * (clinical_risk_group - 1)
  )

  event_time <- rexp(n, rate = event_rate)
  censor_time <- runif(n, min = 180, max = 1200)

  data.frame(
    follow_up_days = pmin(event_time, censor_time),
    event_status = as.integer(event_time <= censor_time),
    biomarker_score = biomarker_score,
    clinical_risk_group = clinical_risk_group
  )
}

discovery_data <- simulate_cohort(500)
validation_data <- simulate_cohort(300)

# For real data, replace the simulation above with, for example:
# discovery_data <- read.csv("data/discovery_data.csv")
# validation_data <- read.csv("data/validation_data.csv")


# ==============================================================================
# 2. Model and timeROC helpers
# ==============================================================================

fit_auc_models <- function(
    data,
    marker_var,
    risk_var,
    time_var,
    event_var,
    evaluation_times
) {
  required_vars <- c(time_var, event_var, marker_var, risk_var)

  analysis_data <- data %>%
    select(all_of(required_vars)) %>%
    filter(complete.cases(.))

  analysis_data[[risk_var]] <- factor(analysis_data[[risk_var]])

  make_cox_model <- function(predictors) {
    coxph(
      reformulate(
        predictors,
        response = paste0("Surv(", time_var, ", ", event_var, ")")
      ),
      data = analysis_data,
      x = TRUE
    )
  }

  make_time_roc <- function(model) {
    timeROC(
      T = analysis_data[[time_var]],
      delta = analysis_data[[event_var]],
      marker = predict(model, newdata = analysis_data, type = "lp"),
      cause = 1,
      times = evaluation_times,
      iid = TRUE
    )
  }

  models <- list(
    biomarker = make_cox_model(marker_var),
    clinical_risk = make_cox_model(risk_var),
    combined = make_cox_model(c(marker_var, risk_var))
  )

  list(
    data = analysis_data,
    models = models,
    rocs = lapply(models, make_time_roc)
  )
}

extract_roc_coordinates <- function(roc_object, time_index, label) {
  data.frame(
    false_positive_rate = roc_object$FP[, time_index],
    true_positive_rate = roc_object$TP[, time_index],
    model = label
  )
}

format_p_value <- function(p) {
  ifelse(
    is.na(p),
    "NA",
    ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
  )
}


# ==============================================================================
# 3. Run discovery and validation analyses
# ==============================================================================

evaluation_times <- c(365, 730)
plot_time <- 365

discovery_results <- fit_auc_models(
  data = discovery_data,
  marker_var = "biomarker_score",
  risk_var = "clinical_risk_group",
  time_var = "follow_up_days",
  event_var = "event_status",
  evaluation_times = evaluation_times
)

validation_results <- fit_auc_models(
  data = validation_data,
  marker_var = "biomarker_score",
  risk_var = "clinical_risk_group",
  time_var = "follow_up_days",
  event_var = "event_status",
  evaluation_times = evaluation_times
)

time_index_discovery <- match(
  plot_time,
  discovery_results$rocs$biomarker$times
)
time_index_validation <- match(
  plot_time,
  validation_results$rocs$biomarker$times
)


# ==============================================================================
# 4. Paired AUC comparisons
# ==============================================================================

# timeROC::compare() requires both prediction markers to be evaluated on the
# same subjects. That condition is satisfied within each cohort because the
# models use the same complete-case analysis dataset.

discovery_comparison <- timeROC::compare(
  discovery_results$rocs$combined,
  discovery_results$rocs$clinical_risk,
  adjusted = FALSE
)

validation_comparison <- timeROC::compare(
  validation_results$rocs$combined,
  validation_results$rocs$clinical_risk,
  adjusted = FALSE
)

discovery_p <- discovery_comparison$p_values_AUC[time_index_discovery]
validation_p <- validation_comparison$p_values_AUC[time_index_validation]

auc_results <- data.frame(
  cohort = c("Discovery", "Validation"),
  evaluation_time = plot_time,
  biomarker_AUC = c(
    discovery_results$rocs$biomarker$AUC[time_index_discovery],
    validation_results$rocs$biomarker$AUC[time_index_validation]
  ),
  clinical_risk_AUC = c(
    discovery_results$rocs$clinical_risk$AUC[time_index_discovery],
    validation_results$rocs$clinical_risk$AUC[time_index_validation]
  ),
  combined_AUC = c(
    discovery_results$rocs$combined$AUC[time_index_discovery],
    validation_results$rocs$combined$AUC[time_index_validation]
  ),
  AUC_difference = c(
    discovery_results$rocs$combined$AUC[time_index_discovery] -
      discovery_results$rocs$clinical_risk$AUC[time_index_discovery],
    validation_results$rocs$combined$AUC[time_index_validation] -
      validation_results$rocs$clinical_risk$AUC[time_index_validation]
  ),
  p_value = c(discovery_p, validation_p),
  p_value_formatted = c(
    format_p_value(discovery_p),
    format_p_value(validation_p)
  )
)

print(auc_results)
write.csv(auc_results, "time_dependent_auc_results.csv", row.names = FALSE)


# ==============================================================================
# 5. Six-curve ROC figure
# ==============================================================================

make_label <- function(cohort, model, auc) {
  sprintf("%s: %s (AUC = %.3f)", cohort, model, auc)
}

labels <- c(
  discovery_biomarker = make_label(
    "Discovery", "Biomarker", auc_results$biomarker_AUC[1]
  ),
  discovery_combined = make_label(
    "Discovery", "Biomarker + clinical risk", auc_results$combined_AUC[1]
  ),
  discovery_risk = make_label(
    "Discovery", "Clinical risk", auc_results$clinical_risk_AUC[1]
  ),
  validation_biomarker = make_label(
    "Validation", "Biomarker", auc_results$biomarker_AUC[2]
  ),
  validation_combined = make_label(
    "Validation", "Biomarker + clinical risk", auc_results$combined_AUC[2]
  ),
  validation_risk = make_label(
    "Validation", "Clinical risk", auc_results$clinical_risk_AUC[2]
  )
)

roc_plot_data <- bind_rows(
  extract_roc_coordinates(
    discovery_results$rocs$biomarker,
    time_index_discovery,
    labels[["discovery_biomarker"]]
  ),
  extract_roc_coordinates(
    discovery_results$rocs$combined,
    time_index_discovery,
    labels[["discovery_combined"]]
  ),
  extract_roc_coordinates(
    discovery_results$rocs$clinical_risk,
    time_index_discovery,
    labels[["discovery_risk"]]
  ),
  extract_roc_coordinates(
    validation_results$rocs$biomarker,
    time_index_validation,
    labels[["validation_biomarker"]]
  ),
  extract_roc_coordinates(
    validation_results$rocs$combined,
    time_index_validation,
    labels[["validation_combined"]]
  ),
  extract_roc_coordinates(
    validation_results$rocs$clinical_risk,
    time_index_validation,
    labels[["validation_risk"]]
  )
)

roc_plot_data$model <- factor(
  roc_plot_data$model,
  levels = unname(labels)
)

roc_plot <- ggplot(
  roc_plot_data,
  aes(false_positive_rate, true_positive_rate, color = model)
) +
  geom_line(linewidth = 1.1) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    color = "grey50"
  ) +
  scale_color_manual(
    values = c(
      "#2166AC", "#004529", "#543005",
      "#B2182B", "darkorange", "#762A83"
    )
  ) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  coord_equal() +
  labs(
    title = paste("Time-dependent ROC at day", plot_time),
    x = "1 - Specificity",
    y = "Sensitivity",
    color = NULL
  ) +
  theme_classic(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

print(roc_plot)

ggsave(
  "time_dependent_roc_example.png",
  roc_plot,
  width = 20,
  height = 18,
  units = "cm",
  dpi = 600,
  bg = "white"
)
