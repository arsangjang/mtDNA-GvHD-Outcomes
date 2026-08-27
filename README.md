Project overview

Graft-versus-host disease (GvHD) remains a major contributor to morbidity and non-relapse mortality after allogeneic hematopoietic cell transplantation (allo-HCT). To investigate whether mitochondrial DNA (mtDNA) variation contributes to GvHD risk, we examined both recipient and donor mtDNA variation in patients with myelodysplastic neoplasms (MDS) undergoing allo-HCT.

In the discovery cohort, whole-genome sequencing was performed on pre-transplant samples from 494 donor–recipient pairs with MDS from the Center for International Blood and Marrow Transplant Research (CIBMTR). Associations of recipient and donor mtDNA variation with acute GvHD grades II–IV (aGvHD24), acute GvHD grades III–IV (aGvHD34), and chronic GvHD (cGvHD) were evaluated using Cox proportional hazards models. To assess the clinical relevance of mtDNA variation for risk prediction, we further compared the performance of models based on the CIBMTR acute GvHD Risk Index with and without mtDNA-derived features, and quantified discrimination using the area under the receiver operating characteristic curve (AUC). Key findings were subsequently evaluated in an independent validation cohort of 295 patients with MDS using targeted deep mtDNA sequencing.

This repository contains the code used for mtDNA processing, gene- and variant-level analyses, association testing, prediction analyses, and figure generation for the study.

Time-dependent AUC example

This repository example demonstrates how to evaluate the incremental predictive
value of a biomarker beyond a clinical risk classification in discovery and
validation cohorts with censored time-to-event outcomes.

The script fits three Cox proportional hazards models in each cohort:

1. Biomarker alone
2. Clinical risk group alone
3. Biomarker plus clinical risk group

It estimates cumulative/dynamic time-dependent ROC curves and AUCs with the
`timeROC` package. Within each cohort, it performs a paired test comparing the
AUC of the clinical-risk-only model with that of the combined model.

## Requirements

Install the required R packages:

```r
install.packages(c("survival", "timeROC", "ggplot2", "dplyr"))
```

## Run

From R or RStudio:

```r
source("time_dependent_auc_example.R")
```

The example is fully reproducible and generates simulated data. To analyze real
data, replace the simulation section with code that imports two de-identified
cohort files containing the following variables:

| Variable | Description |
|---|---|
| `follow_up_days` | Observed event or censoring time |
| `event_status` | Event indicator (`1` = event, `0` = censored) |
| `biomarker_score` | Continuous or ordinal biomarker |
| `clinical_risk_group` | Clinical risk category |

## Outputs

- `time_dependent_auc_results.csv`: cohort-specific AUCs, AUC differences, and p-values
- `time_dependent_roc_example.png`: six-curve discovery/validation ROC figure

## Statistical comparison

The paired comparison is:

```text
Clinical risk alone versus biomarker + clinical risk
