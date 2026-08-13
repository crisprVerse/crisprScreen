#' @title Generate plot of barcode distributions on a log2 scale (static plot)
#' @description Generate plot of barcode distributions on a log2 scale (static plot).
#' @param object Either a numeric matric or a SummarizedExperiment
#' @param log.scale Should the data be log-transformed? TRUE by default.
#' @param sort Should samples be drawn by color order? 
#' @param col Optional character vector specifying sample color. Must be 
#' @param col.group.field Optional string indicating which column of colData(object) should be
#' used to color densities. 
#' @param lty Optional character vector specifying sample line type. Must be 
#' of the same length as the number of columns in the provided object.
#' @param xlim Usual xlim argument in base plot
#' @param ylim Usual ylim argument in base plot
#' @param xlab Usual xlab argument in base plot
#' @param ylab Usual ylab argument in base plot
#' @param main Usual main argument in base plot
#' @param legend.where Where should be the legend?
#' @param assay_index Integer specifying which assay should be used. 1 by default. 
#' @param ... Other arguments to passe to base plot function
#' @return Returns nothing. A plot is produced as a side effect. 
#' @examples 
#' \dontrun{
#'  data("seExample")
#'  drawDensities(seExample)
#' }
#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SummarizedExperiment assays
#' @importFrom matrixStats rowMedians
#' @importFrom stats density
#' @importFrom graphics lines plot legend
drawDensities <- function(object,
    log.scale=TRUE, 
    sort=TRUE, 
    col=NULL,
    col.group.field="Group",
    lty=NULL,
    xlim=NULL, 
    ylim=NULL, 
    xlab="Log2 sgRNA Counts", 
    ylab="Density", 
    main="",
    legend.where="topleft",
    assay_index=1,
    ...
){
    if (is.null(lty)){
    	lty <- rep(1, ncol(object))
    }
    n <- ncol(object)
	#if (class(object)=="ExpressionSet") counts <- exprs(object)
    if (is(object, "SummarizedExperiment")){
        xs <- assays(object)
        if (assay_index>length(xs)){
            stop("Assay index is greated than the number of assays. ")
        } else {
            counts <- as.matrix(xs[[assay_index]])
        }
    } else {
        counts <- object
    }
	if (log.scale){
        counts <- log2(counts+1)
    } 
    if (is.null(col)){
        pheno <- colData(object)
        if (!col.group.field %in% colnames(pheno)){
            stop("col.group.field is not found in colData(object)")
        }
        key     <- data.frame(label=pheno[[col.group.field]])
        key$col <- as.numeric(as.factor(key$label))
        col <- key$col
    } 
	o <- rep(1, n)
	if (sort){
        o <- order(col)    
    }
    densities <- lapply(1:n, function(i){
        density(counts[,o][,i], na.rm=TRUE)
    })

    # Setting up axes limits:
    minx <- min(unlist(lapply(densities, function(x) min(x$x))))
    maxx <- max(unlist(lapply(densities, function(x) max(x$x))))
    miny <- min(unlist(lapply(densities, function(x) min(x$y))))
    maxy <- max(unlist(lapply(densities, function(x) max(x$y))))
    rangex <- abs(maxx-minx)
    rangey <- abs(maxy-miny)
    epsx <- 0.1*rangex
    epsy <- 0.1*rangey

    if (is.null(xlim)){
        xlim=c(minx-epsx,maxx+epsx)
    } 
    if (is.null(ylim)){
      ylim=c(miny,maxy+epsy)  
    } 

	plot(densities[[1]], col="white", bty="L", 
        xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, main=main,...)
	for (i in 1:ncol(counts)){
		lines(densities[[i]], col=col[o][i], lty=lty[o][i],...)
	}
    legend(legend.where, col=unique(key$col), legend=unique(key$label), lty=1, bty="n")
}













