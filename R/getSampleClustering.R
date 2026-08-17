#' @title Cluster samples using log-ratio profiles
#'
#' @description
#' Performs hierarchical clustering of samples based on feature-level log
#' ratios relative to a reference condition. The function first computes log
#' ratios using \code{\link{getLogRatios}}, removes samples belonging to the
#' reference condition, selects the most variable features across the remaining
#' samples, and calculates a sample-by-sample similarity matrix. Hierarchical
#' clustering is then performed using a distance defined as
#' \code{1 - similarity}.
#'
#' @param se A \code{SummarizedExperiment} object containing feature-level
#'   measurements.
#' @param condition.field Character string specifying the column in
#'   \code{colData(se)} containing the experimental condition.
#' @param reference.level Character string specifying the reference condition
#'   used to calculate log ratios.
#' @param replicate.field Character string specifying the column in
#'   \code{colData(se)} identifying biological or technical replicates used
#'   when calculating log ratios.
#' @param how.many Integer specifying the maximum number of most variable
#'   features to use for sample clustering. Default is \code{5000}. If fewer
#'   eligible features are available, all eligible features are used.
#' @param method Character string specifying the correlation method used to
#'   calculate sample similarity. Passed to \code{\link[stats]{cor}}.
#'   Common choices are \code{"pearson"}, \code{"spearman"}, and
#'   \code{"kendall"}. Default is \code{"pearson"}.
#'
#' @details
#' Log ratios are first calculated using \code{\link{getLogRatios}}. Samples
#' belonging to the reference condition are subsequently excluded from the
#' clustering analysis.
#'
#' Features containing non-finite values or having zero variance across samples
#' are removed. Among the remaining features, the \code{how.many} features with
#' the largest variance are retained.
#'
#' Sample similarity is calculated by correlating the selected feature profiles.
#' The similarity matrix is converted to a distance matrix as
#' \code{1 - similarity}, and hierarchical clustering is performed using
#' average linkage via \code{\link[stats]{hclust}}.
#'
#' @return A named list containing:
#' \describe{
#'   \item{\code{similarity}}{A sample-by-sample similarity matrix.}
#'   \item{\code{distance}}{A \code{dist} object containing distances defined
#'     as \code{1 - similarity}.}
#'   \item{\code{clustering}}{An \code{hclust} object containing the
#'     hierarchical clustering result.}
#'   \item{\code{features}}{Character vector containing the feature identifiers
#'     used for clustering.}
#' }
#'
#' @seealso
#' \code{\link{getLogRatios}}, \code{\link[stats]{cor}},
#' \code{\link[stats]{hclust}}
#'
#' @importFrom matrixStats rowVars
#' @importFrom stats as.dist hclust
#'
#' @export
getSampleClustering <- function(se,
                                condition.field,
                                reference.level,
                                replicate.field,
                                how.many = 5000,
                                method = "pearson") {

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

    ## Calculate feature variance across samples
    vars <- rowVars(Y)

    ## Remove invariant features
    keep <- is.finite(vars) & vars > 0
    Y <- Y[keep, , drop = FALSE]
    vars <- vars[keep]

    ## Select most variable features
    how.many <- min(how.many, nrow(Y))
    wh <- order(vars, decreasing = TRUE)[seq_len(how.many)]
    Y <- Y[wh, , drop = FALSE]

    ## Sample-by-sample correlation
    similarity <- stats::cor(Y, method = method)

    ## Convert correlation into a distance
    distance <- as.dist(1 - similarity)

    ## Hierarchical clustering
    clustering <- hclust(
        distance,
        method = "average"
    )

    return(list(
        similarity = similarity,
        distance = distance,
        clustering = clustering,
        features = rownames(Y)
    ))
}
