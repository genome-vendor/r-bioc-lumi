bgAdjust <- function(lumiBatch, probs=0.5, ...) {
	if (!is(lumiBatch, 'LumiBatch')) stop('The object should be class "LumiBatch"!')
	## by subtract an offset, which is estimated based on the quantile of the control probes
	if (min(exprs(lumiBatch)) <= 1) {
		print('The data has already been background adjusted!')
		return(lumiBatch)
	}
	control <- lumiBatch@controlData
	if (is.null(control) || nrow(control) == 0) {
		print('There is no control probe information in the LumiBatch object!\n No background adjustment will be performed.')
		return(lumiBatch)
	}
	control <- control[, sampleNames(lumiBatch)]
	probeType <- rownames(control)
	if ('negative' %in% probeType) control <- control[probeType == 'negative',]
	quantile.ctrl <- apply(control, 2, quantile, probs=probs, ...)
	exprs(lumiBatch) <- exprs(lumiBatch) - matrix(rep(1, nrow(lumiBatch)), ncol=1) %*% quantile.ctrl
	return(lumiBatch)
}