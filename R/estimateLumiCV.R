`estimateLumiCV` <-
function(x.lumi, type=c('measurement', 'probe'), ifPlot=FALSE, ...) {

	type <- match.arg(type)
	if (type == 'measurement') {
		if (!is(x.lumi, 'LumiBatch')) stop('The object should be class "LumiBatch"!')
		std <- se.exprs(x.lumi)
		if (is.null(std)) stop('The object does not include "se.exprs" matrix in assayData slot!')
		cv <- std / exprs(x.lumi)
		rownames(cv) <- rownames(std)
		colnames(cv) <- colnames(std)
		if (ifPlot) {
			plotDensity(cv, xlab='coefficient of variance', ...)
			return(invisible(cv))	
		}
	} else {
		if (is(x.lumi, 'ExpressionSet')) x.lumi <- exprs(x.lumi)
		cv <- apply(x.lumi, 1, function(x) sd(x)/mean(x))
		names(cv) <- rownames(x.lumi)
	}
	return(cv)
}

