# m6APrediction

## Overview

`m6APrediction` is an R package for predicting N6-methyladenosine (m6A) sites from sequence and genomic features.

The package uses a fitted random forest model to calculate the predicted probability of an m6A site and classify each site as either `Positive` or `Negative`.

The package provides two exported functions:

- `prediction_multiple()` predicts m6A status for multiple sites.
- `prediction_single()` predicts m6A status for one site.

## Installation

You can install the development version from GitHub using:

```r
install.packages("remotes")
remotes::install_github("JiahanZhang-xjtlu/m6APrediction")
```

## Usage

Predict multiple sites

```r
library(m6APrediction)

ml_fit <- readRDS(system.file(
  "extdata",
  "rf_fit.rds",
  package = "m6APrediction",
  mustWork = TRUE
))

example_df <- read.csv(system.file(
  "extdata",
  "m6A_input_example.csv",
  package = "m6APrediction",
  mustWork = TRUE
))

batch_result <- prediction_multiple(
  ml_fit = ml_fit,
  feature_df = example_df,
  positive_threshold = 0.6
)

head(batch_result)
```

The returned data frame contains the original input features and two
additional columns:

- predicted_m6A_prob: predicted probability of an m6A site.
- predicted_m6A_status: classification as Positive or Negative.

## Predict a single site

```r
library(m6APrediction)

ml_fit <- readRDS(system.file(
  "extdata",
  "rf_fit.rds",
  package = "m6APrediction",
  mustWork = TRUE
))

single_result <- prediction_single(
  ml_fit = ml_fit,
  gc_content = 0.6,
  RNA_type = "mRNA",
  RNA_region = "CDS",
  exon_length = 12,
  distance_to_junction = 5,
  evolutionary_conservation = 0.8,
  DNA_5mer = "GGACA",
  positive_threshold = 0.5
)

single_result
```

The result is a named character vector containing:

- predicted_m6A_prob
- predicted_m6A_status

## Input features

The prediction functions use the following features:

- GC content
- RNA type
- RNA region
- exon length
- distance to the nearest splice junction
- evolutionary conservation
- five-base DNA sequence context

The exon length and distance-to-junction variables should use the same
log2-transformed format as the model training data.

## Model performance

The following figure shows the model performance on the evaluation data. The left panel shows the receiver operating characteristic (ROC) curve, and the right panel shows the precision-recall curve. The evaluation data contained 2,273 positive observations and 7,338
negative observations.

## Model resources

The package includes the following resources:

inst/extdata/rf_fit.rds
inst/extdata/m6A_input_example.csv

The example CSV file contains 100 example rows for demonstrating package usage. It is not an independent test set for evaluating model accuracy.

## License
This project is licensed under the MIT License.