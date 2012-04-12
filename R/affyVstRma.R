affyVstRma <- function (afbatch, bgcorrect.method = 'none', bgcorrect.param = list(), VST.param = list(), verbose = TRUE, ...) 
{
	if (is.null(bgcorrect.method)) bgcorrect.method <- 'none'
	if (is.character(afbatch)) {
		if (!all(file.exists(afbatch))) stop('The input is not a valid file list or an AffyBatch object!')
		fileList <- afbatch
		afbatch <- NULL
		for (i in 1:length(fileList)) {
			file.i <- fileList[i]
			cat(paste('Processing', file.i, '...\n'))
			afbatch.i <- ReadAffy(filenames=file.i, sd=TRUE)

			if (bgcorrect.method != 'none') {
				if (verbose)  cat("background correcting...")
				afbatch.i <- do.call("bg.correct", c(alist(afbatch.i, method = bgcorrect.method), bgcorrect.param))
				if (verbose)  cat("done.\n")
			}
			if (verbose) cat("Variance stabilizing ...\n")
	# browser()
			afbatch.i <- do.call("lumiT", c(alist(afbatch.i), VST.param))

			if (verbose) cat("done.\n")
			if (i == 1) {
				afbatch <- afbatch.i
			} else {
				afbatch <- merge.AffyBatch(afbatch, afbatch.i)
			}
		}
	} else {
		if (!is(afbatch, 'AffyBatch')) stop('The input is not an AffyBatch object or a valid file list!')
		if (bgcorrect.method != 'none') {
			if (verbose)  cat("background correcting...")
			afbatch <- do.call("bg.correct", c(alist(afbatch, method = bgcorrect.method), bgcorrect.param))
			if (verbose)  cat("done.\n")
		}
		if (verbose) cat("Variance stabilizing ...\n")
		afbatch <- do.call("lumiT", c(alist(afbatch), VST.param))
		if (verbose) cat("done.\n")
	}
	## To avoid double transformation in the rma function
	exprs(afbatch) <- 2^exprs(afbatch)

	if (verbose) cat("RMA normalization ...\n")
	eset <- rma(afbatch, background=FALSE, ...)
	if (verbose) cat("done.\n")
	return(eset)
}
