# m6APrediction

## Overview

`m6APrediction` is an R package for predicting N6-methyladenosine
(m6A) sites from sequence and genomic features.

The package uses a fitted random forest model to calculate the predicted
probability of an m6A site and classify each site as either `Positive`
or `Negative`.

The package provides two exported functions:

- `prediction_multiple()` predicts m6A status for multiple sites.
- `prediction_single()` predicts m6A status for one site.

## Installation

You can install the development version from GitHub using:

```r
install.packages("remotes")
remotes::install_github("JiahanZhang-xjtlu/m6APrediction")
