#' @title Get ROC curve for classifying essential genes
#' @description Get receiver operating characteristic (ROC) curve for 
#' classifying essential genes.
#' @param x A numeric vector
#' @param genes A character vector specifying genes associated with x
#' @param rev Should negative values represent signal? TRUE by default. 
#' @param essentialGenes Character vector of essential gene symbols.
#' @param nonEssentialGenes Character vector of non-essential gene symbols.
#' @param pseudo Should pseudo ROC curves be computed instead?
#'     FALSE by default.
#' @return A data.frame with 2 columns:
#'     \code{specificity} and \code{sensitivity}.
#' @export
#' @importFrom stats approxfun
getEssentialityROC <- function(x,
                               genes,
                               rev=TRUE, 
                               essentialGenes,
                               nonEssentialGenes,
                               pseudo=FALSE
){
    if (rev){
        x <- -x
    }
    if (length(x)!=length(genes)){
        stop("x must have the same length as genes")
    } 
  

    isEssential <- genes %in% essentialGenes
    # In the case genes are mouse genes:
    if(sum(isEssential)==0){
        stop("[getEssentialityROC] Cannot find essential ",
             "genes among provided genes.")
    }
    if (!pseudo){
        isNonEssential <- genes %in% nonEssentialGenes
        if (sum(isNonEssential)==0){
            stop("[getEssentialityROC] Cannot find non-essential genes",
                 " among provided genes. Try with pseudo=TRUE.")
        }
    }
    if (pseudo){
        tmp <- .roc_fast(x, as.numeric(isEssential))
    } else {
        wh <- isEssential | isNonEssential
        tmp <- .roc_fast(x[wh], as.numeric(isEssential[wh]))
    }
    out <- data.frame(specificity=tmp$specificities,
                      sensitivity=tmp$sensitivities)
    return(out)
}


.roc_fast <- function(probs, class){
    # Removing non-finite values:
    wh <- which(is.finite(probs))
    if (length(wh)==0){
        stop("[getEssentialityROC] No finite values provided.")
    }
    probs <- probs[wh]
    class <- class[wh]
    # Calculating ROC:
    class_sorted <- class[order(probs, decreasing=T)]
    TPR <- cumsum(class_sorted) / sum(class)
    FPR <- cumsum(class_sorted == 0) / sum(class == 0)
    return(list(sensitivities=TPR,
                specificities=1-FPR))
}




