#' @title Get principal components
#'
#' @description
#' Performs principal component analysis on log-ratio profiles derived from a
#' SummarizedExperiment object. Log-ratios are calculated relative to a
#' reference condition using \code{getLogRatios()}. PCA is performed across
#' samples using the most variable features.
#'
#' @param se A SummarizedExperiment object.
#' @param condition.field Character string specifying the column in
#'   \code{colData(se)} containing the experimental condition.
#' @param reference.level Character string specifying the reference condition
#'   used to calculate log-ratios.
#' @param replicate.field Character string specifying the replicate column in
#'   \code{colData(se)}.
#' @param how.many Maximum number of most-variable features to use for PCA.
#'   Default is 5000.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{pc1}: sample scores for principal component 1.
#'   \item \code{pc2}: sample scores for principal component 2.
#'   \item \code{perc1}: percentage of variance explained by PC1.
#'   \item \code{perc2}: percentage of variance explained by PC2.
#'   \item \code{scores}: matrix containing all sample PC scores.
#'   \item \code{variance.explained}: percentage of variance explained by
#'     each principal component.
#'   \item \code{pca}: the complete object returned by \code{prcomp()}.
#' }
#'
#' @importFrom stats prcomp
#' @importFrom matrixStats rowVars
#' @export
getPCs <- function(se,
                   condition.field,
                   reference.level,
                   replicate.field,
                   how.many = 5000) {

    seRatios <- getLogRatios(
        se,
        condition.field = condition.field,
        reference.level = reference.level,
        replicate.field = replicate.field
    )
    pheno <- colData(seRatios)
    seRatios <- seRatios[, pheno[[condition.field]]!=reference.level]
    Y <- assays(seRatios)[[1]]

    ## Remove features containing non-finite values
    keep <- rowSums(!is.finite(Y)) == 0
    Y <- Y[keep, , drop = FALSE]

    ## Calculate variance of each feature across samples
    vars <- rowVars(Y)

    ## Remove features with no variation
    keep <- is.finite(vars) & vars > 0
    Y <- Y[keep, , drop = FALSE]
    vars <- vars[keep]

    ## Select the most variable features
    how.many <- min(how.many, nrow(Y))
    wh <- order(vars, decreasing = TRUE)[seq_len(how.many)]
    Y <- Y[wh, , drop = FALSE]

    ## PCA:
    ## rows = samples
    ## columns = CRISPR features
    pca <- prcomp(
        t(Y),
        center = TRUE,
        scale. = FALSE
    )

    ## Percentage variance explained
    variance.explained <- 100 * pca$sdev^2 / sum(pca$sdev^2)

    out <- list(
        pc1 = pca$x[, 1],
        pc2 = pca$x[, 2],
        perc1 = round(variance.explained[1], 1),
        perc2 = round(variance.explained[2], 1),
        scores = pca$x,
        variance.explained = variance.explained,
        pca = pca,
        samples=colnames(seRatios)
    )
    return(out)
}



