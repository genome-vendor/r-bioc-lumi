`rsn` <-
function(x.lumi, targetArray=NULL, excludeFold=2, span=0.03, ifPlot=FALSE,...) {
	if (is(x.lumi, 'ExpressionSet')) {
		# x.lumi is a lumi object
		exprs <- exprs(x.lumi)
	} else if (is.numeric(x.lumi)) {
		exprs <- as.matrix(x.lumi)
	} else {
		stop('The object should be a matrix or class "ExpressionSet" inherited!')
	}

	externalTarget <- FALSE
	if (!is.null(targetArray)) {
		## check the format of the targetArray
		if (is(targetArray, 'ExpressionSet')) {
			targetArray <- exprs(targetArray)[,1]
		} 
		if (length(targetArray) > 1) {
			if (length(targetArray) != nrow(exprs)) stop('targetArray should be an index or a vector has the same length as other samples.')
			exprs <- cbind(targetArray, exprs)
			targetArray <- 1
			externalTarget <- TRUE
		}
		if (is.numeric(targetArray)) {
			if (targetArray < 1 || targetArray > nrow(exprs)) {
				warning('The provided targetArray is invalid and will be set as NULL!')
				targetArray <- NULL
			}
		} else if (is.character(targetArray)) {
			if (!(targetArray %in% colnames(exprs))) {
				warning('The provided targetArray is invalid and will be set as NULL!')
				targetArray <- NULL
			}
		} else {
			warning('The provided targetArray is invalid and will be set as NULL!')
			targetArray <- NULL
		}
	}

	## check whether the data was variance stabilized.
	if (max(exprs, na.rm=TRUE) > 100) {
		log2Trans <- FALSE
		
		offset <- ifelse(min(exprs) < 0, min(exprs) - 1, 0)
		exprs <- log2(exprs - offset)
	} else {
		log2Trans <- TRUE
	}

	## Define interal function 
	pairwiseN <- function(ind, exprs, targetArray, exprs0=NULL, ifPlot=FALSE) {
		cat(as.character(Sys.time()), ", processing array ", ind, "\n")

		# normal array ind against targetArray      
		u1 <- exprs[,ind]
		u2 <- exprs[, targetArray]
		u1.original <- u1
		u2.original <- u2

		## define window functions
		win <- function(x, sigma=1, type=c('gaussian', 'tricube', 'bisqure')) {
			type <- match.arg(type)
			ww <- switch(type,
				gaussian = exp(-x^2/(2*sigma^2)) / (sqrt(2*pi)*sigma),
				tricube = (1 - (abs(x)/(3 * sigma))^3)^3,
				bisquare = (1 - (x/(3 * sigma))^2)^2 )
			return(ww)
		}

		if (ind != targetArray[1] || length(targetArray) > 1) {
			## calculate the weights based on the fold change 
			if (!is.null(exprs0)) {
				fd <- exprs0[,ind] - exprs0[, targetArray]
				if (!is.null(excludeFold)) {
					sigma <- log2(excludeFold)/3
				} else {
					sigma <- sd(fd)
				}
 				wt <- win(abs(fd), sigma)
				wt <- wt/max(wt)
			} else {
				wt <- rep(length(u1))
			}

			u1.normalized <- monoSmu(u1, u2, newX=u1.original, span=span, ifPlot=ifPlot, wt=wt, rotate=TRUE,
					xlab=paste("array", ind), ylab=paste("array", targetArray), ...)
		} else {
			u1.normalized <- u1
		}
		return(u1.normalized)
	}

	if (ifPlot) par(mfrow=c(2,2))

	## do quantile normalization for the purpose of estimating fold change
	## Based on the estimated fold change, we can down-weight the differentiated genes.
	if (!is.null(excludeFold)) {
		exprs0 <- normalize.quantiles(exprs)
	} else {
		exprs0 <- NULL
	}

	if (is.null(targetArray)) {
		# find the sample which is the most similar to the mean profile of all samples,
		meanProfile <- apply(exprs, 1, mean)
		targetArray <- which.min(abs(colSums(exprs - meanProfile)))
	}

	nArray <- ncol(exprs)
	normalized <- lapply(1:nArray, FUN=pairwiseN, exprs=exprs,
					exprs0=exprs0, targetArray=targetArray, ifPlot=ifPlot)
	normalized <- matrix(unlist(normalized), ncol=nArray, byrow=FALSE)	
	colnames(normalized) <- colnames(exprs)
	rownames(normalized) <- rownames(exprs)

	## if the targetArray is an external vector, it will be removed from the normalized data.
	if (externalTarget) normalized <- normalized[,-1]
	## transformed as original scale in not log2transformed
	if (!log2Trans) normalized <- 2^normalized
	
	if (is(x.lumi, 'ExpressionSet')) {
		exprs(x.lumi) <- normalized
		if (is.numeric(targetArray)) {
			if (!is.null(attr(x.lumi, 'vstParameter'))) {
				attr(x.lumi, 'vstParameter') <- attr(x.lumi, 'vstParameter')[targetArray,]
				attr(x.lumi, 'transformFun') <- attr(x.lumi, 'transformFun')[targetArray]
				attr(x.lumi, 'targetArray') <- sampleNames(x.lumi)[targetArray]
			}
		} else {
			if (!is.null(attr(targetArray, 'vstParameter'))) {
				attr(x.lumi, 'vstParameter') <- attr(targetArray, 'vstParameter')
				attr(x.lumi, 'transformFun') <- attr(targetArray, 'transformFun')
			}
			if (is(targetArray, 'ExpressionSet'))
				attr(x.lumi, 'targetArray') <- sampleNames(targetArray)
		}
	} else {
		x.lumi <- normalized
	}
	
	return(x.lumi)
}
