#' Encode DNA sequences into nucleotide factors
#'
#' Converts DNA sequences into a data frame containing one factor column
#' for each nucleotide position.
#'
#' @param dna_strings Character vector of DNA sequences with equal length.
#' @return A data frame containing nucleotide positions encoded as factors.
#' @importFrom stats predict
#'
dna_encoding <- function(dna_strings){
  nn <- nchar(dna_strings[1])
  seq_m <- matrix(unlist(strsplit(dna_strings, "")), ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  return(seq_df)
}
#' Predict m6A status for multiple sites
#'
#' Uses a fitted random forest model to predict m6A probabilities
#' and classification status for multiple input sites.
#'
#' @param ml_fit A fitted random forest model.
#' @param feature_df A data frame containing the required input features.
#' @param positive_threshold Numeric classification threshold.
#'   Scores strictly greater than this value are classified as Positive.
#' @return The input data frame with predicted_m6A_prob and
#'   predicted_m6A_status columns added.
#' @import randomForest
#' @importFrom stats predict
#' @export
#' @examples
#' ml_fit <- readRDS(system.file(
#'   "extdata",
#'   "rf_fit.rds",
#'   package = "m6APrediction",
#'   mustWork = TRUE
#' ))
#'
#' example_df <- read.csv(system.file(
#'   "extdata",
#'   "m6A_input_example.csv",
#'   package = "m6APrediction",
#'   mustWork = TRUE
#' ))
#'
#' example_df <- example_df[1:5, ]
#'
#' prediction_multiple(
#'   ml_fit = ml_fit,
#'   feature_df = example_df,
#'   positive_threshold = 0.6
#' )
#'
prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5){
  stopifnot(all(c("gc_content", "RNA_type", "RNA_region", "exon_length", "distance_to_junction",
                  "evolutionary_conservation", "DNA_5mer") %in% colnames(feature_df))) # Check errors if incorrect column names of input data.frame
  feature_df$RNA_type <- factor(feature_df$RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene"))
  feature_df$RNA_region <- factor(feature_df$RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR"))
  infered_prob <- predict(ml_fit, newdata = cbind(feature_df, dna_encoding(feature_df$DNA_5mer)), type = "prob")
  feature_df$predicted_m6A_prob <- infered_prob[,"Positive"]
  feature_df$predicted_m6A_status <- "Negative"
  feature_df$predicted_m6A_status[feature_df$predicted_m6A_prob > positive_threshold] <- "Positive"
  return(feature_df) # Return a data.frame with predicted m6A prob and predicted m6A status
}
#' Predict m6A status for a single site
#'
#' Predicts the m6A probability and classification status for one site.
#'
#' @param ml_fit A fitted random forest model.
#' @param gc_content Numeric GC-content value.
#' @param RNA_type Character RNA type.
#' @param RNA_region Character RNA region.
#' @param exon_length Numeric log2-transformed exon length.
#' @param distance_to_junction Numeric log2-transformed distance to
#'   the nearest junction.
#' @param evolutionary_conservation Numeric conservation score.
#' @param DNA_5mer Character five-base sequence with A in the centre.
#' @param positive_threshold Numeric classification threshold.
#'   Scores strictly greater than this value are classified as Positive.
#' @return A named character vector containing the predicted probability
#'   and predicted status.
#' @import randomForest
#' @importFrom stats predict
#' @export
#' @examples
#' ml_fit <- readRDS(system.file(
#'   "extdata",
#'   "rf_fit.rds",
#'   package = "m6APrediction",
#'   mustWork = TRUE
#' ))
#'
#' prediction_single(
#'   ml_fit = ml_fit,
#'   gc_content = 0.6,
#'   RNA_type = "mRNA",
#'   RNA_region = "CDS",
#'   exon_length = 12,
#'   distance_to_junction = 5,
#'   evolutionary_conservation = 0.8,
#'   DNA_5mer = "GGACA",
#'   positive_threshold = 0.5
#' )
#'
prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region, exon_length,
                              distance_to_junction, evolutionary_conservation, DNA_5mer, positive_threshold = 0.5){
  pred_outcome <- prediction_multiple(ml_fit = ml_fit,
                                      data.frame(gc_content = gc_content,
                                                 RNA_type = factor(RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene")),
                                                 RNA_region = factor(RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR")),
                                                 exon_length = exon_length,
                                                 distance_to_junction = distance_to_junction,
                                                 evolutionary_conservation = evolutionary_conservation,
                                                 DNA_5mer = DNA_5mer),
                                      positive_threshold = positive_threshold)
  returned_vector <- c(pred_outcome$predicted_m6A_prob, pred_outcome$predicted_m6A_status)
  names(returned_vector) <- c("predicted_m6A_prob", "predicted_m6A_status")
  return(returned_vector) # Return a named vector with predicted m6A prob and predicted m6A status
}
