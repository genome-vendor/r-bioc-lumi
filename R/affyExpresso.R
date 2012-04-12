affyExpresso <- function (afbatch, bg.correct = TRUE, bgcorrect.method = NULL, bgcorrect.param = list(), variance.stabilize = TRUE,
    varianceStabilize.method = c("vst", "log2", "cubicRoot"), varianceStabilize.param = list(), normalize = TRUE, normalize.method = NULL, 
    normalize.param = list(), pmcorrect.method = NULL, pmcorrect.param = list(), 
    summary.method = NULL, summary.param = list(), summary.subset = NULL, verbose = TRUE) 
{
    varianceStabilize.method <- match.arg(varianceStabilize.method)
	setCorrections <- function() {
        bioc.opt <- getOption("BioC")
        if (bg.correct) {
            if (is.null(bgcorrect.method)) {
                BGMethods <- bgcorrect.methods
            } else {
                BGMethods <- bgcorrect.method
            }
        } else {
            BGMethods <- "None"
        }
        if (normalize) {
            if (is.null(normalize.method)) {
                normMethods <- normalize.methods(afbatch)
            } else {
                normMethods <- normalize.method
            }
        } else {
            normMethods <- "None"
        }
        if (is.null(pmcorrect.method)) {
            PMMethods <- pmcorrect.methods
        } else {
            PMMethods <- pmcorrect.method
        }
        if (is.null(summary.method)) {
            expMethods <- generateExprSet.methods
        } else {
            expMethods <- summary.method
        }
        corrections <- expressoWidget(BGMethods, normMethods, 
            PMMethods, expMethods, bioc.opt$affy$bgcorrect.method, 
            bioc.opt$affy$normalize.method, bioc.opt$affy$pmcorrect.method, 
            bioc.opt$affy$summary.method)
        if (!is.null(corrections)) {
            if (corrections[["BG"]] != "None") {
                bgcorrect.method <<- corrections[["BG"]]
            }
            if (corrections[["NORM"]] != "None") {
                normalize.method <<- corrections[["NORM"]]
            }
            if (corrections[["PM"]] != "None") {
                pmcorrect.method <<- corrections[["PM"]]
            }
            if (corrections[["EXP"]] != "None") {
                summary.method <<- corrections[["EXP"]]
            }
        } else {
            stop("Aborted by user")
        }
    }

    nchips <- length(afbatch)
    if (verbose) {
        if (bg.correct) {
            cat("background correction:", bgcorrect.method, "\n")
        }
		if (variance.stabilize) {
			vstMethod <- ifelse(is.null(varianceStabilize.param$method), 'vst', varianceStabilize.param$method)
			cat("Variance Stabilizing Transform:", vstMethod, "\n")
		}
        if (normalize) {
            cat("normalization:", normalize.method, "\n")
        }
        cat("PM/MM correction :", pmcorrect.method, "\n")
        cat("expression values:", summary.method, "\n")
    }

	if (is.character(afbatch)) {
		if (!all(file.exists(afbatch))) stop('The input is not a valid file list or an AffyBatch object!')
		fileList <- afbatch
		afbatch <- NULL
		for (i in 1:length(fileList)) {
			file.i <- fileList[i]
			cat(paste('Processing', file.i, '...\n'))
			afbatch.i <- ReadAffy(filenames=file.i, sd=TRUE)

		    if (bg.correct) {
		        if (verbose) 
		            cat("background correcting...")
					afbatch.i <- do.call("bg.correct", c(alist(afbatch.i, method = bgcorrect.method), bgcorrect.param))
		        if (verbose) 
		            cat("done.\n")
		    }
			if (variance.stabilize) {
				if (verbose) cat("Variance stabilizing ...\n")
				afbatch.i <- do.call("lumiT", c(alist(afbatch.i, method = varianceStabilize.method), varianceStabilize.param))
				if (verbose) cat("done.\n")
			}
			if (i == 1) {
				afbatch <- afbatch.i
			} else {
				afbatch <- merge.AffyBatch(afbatch, afbatch.i)
			}
		}
	} else {
	    if (bg.correct) {
	        if (verbose) 
	            cat("background correcting...")
	        afbatch <- do.call("bg.correct", c(alist(afbatch, method = bgcorrect.method), 
	            bgcorrect.param))
	        if (verbose) 
	            cat("done.\n")
	    }
		if (variance.stabilize) {
			if (verbose) cat("Variance stabilizing ...\n")
			afbatch <- do.call("lumiT", c(alist(afbatch), varianceStabilize.param))
			if (verbose) cat("done.\n")
		}
	}

    if (normalize) {
        if (verbose) 
            cat("normalizing...")
        afbatch <- do.call("normalize", c(alist(afbatch, normalize.method), 
            normalize.param))
        if (verbose) 
            cat("done.\n")
    }
	## To avoid double transformation in the rma function
	if (summary.method == 'medianpolish' && variance.stabilize)  exprs(afbatch) <- 2^exprs(afbatch)
    eset <- computeExprSet(afbatch, summary.method = summary.method, 
        pmcorrect.method = pmcorrect.method, ids = summary.subset, 
        summary.param = summary.param, pmcorrect.param = pmcorrect.param)
    return(eset)
}