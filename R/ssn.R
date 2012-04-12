`ssn` <-
function(x.lumi, targetArray=NULL, scaling=TRUE, bgMethod=c('density', 'mean', 'median', 'none'), fgMethod=c('mean', 'density', 'median'), ...) {
	if (is(x.lumi, 'ExpressionSet')) {
		expr  <- exprs(x.lumi)
	} else if (is.numeric(x.lumi)) {
		expr <- as.matrix(x.lumi)
	} else {
		stop('The object should be a matrix or class "ExpressionSet" inherited!')
	}
	bgMethod <- match.arg(bgMethod)
	fgMethod <- match.arg(fgMethod)
	if (bgMethod != 'density' && bgMethod == fgMethod) scaling <- FALSE
	if (bgMethod != 'density' && fgMethod == 'density') {
		warning('fgMethod "density" can only be used together with the bgMethod "density"!\n bgMethod is set as "density".')
		bgMethod <- 'density'
	}
	externalTarget <- FALSE
	if (!is.null(targetArray)) {
		## check the format of the targetArray
		if (is(targetArray, 'ExpressionSet')) {
			targetArray <- exprs(targetArray)[,1]
		} 
		if (length(targetArray) > 1) {
			if (length(targetArray) != nrow(expr)) stop('targetArray should be an index or a vector has the same length as other samples.')
			expr <- cbind(targetArray, expr)
			targetArray <- 1
			externalTarget <- TRUE
		}
		if (is.numeric(targetArray)) {
			if (targetArray < 1 || targetArray > nrow(expr)) {
				warning('The provided targetArray is invalid and will be set as NULL!')
				targetArray <- NULL
			}
		} else if (is.character(targetArray)) {
			if (!(targetArray %in% colnames(expr))) {
				warning('The provided targetArray is invalid and will be set as NULL!')
				targetArray <- NULL
			}
		} else {
			warning('The provided targetArray is invalid and will be set as NULL!')
			targetArray <- NULL
		}
	}

	## check whether the data was variance stabilized.
	if (max(expr, na.rm=TRUE) > 100) {
		log2Trans <- FALSE
	} else {
		log2Trans <- TRUE
		expr <- 2^(expr)
	}

	## Estimate the background level and the AUC (mean) after removing the background level
	twoPoint <- apply(expr, 2, function(xx) {
		hh <- hist(xx, 1000, plot=FALSE)
		Th <- hh$breaks[which.max(hh$counts) + 1] * 2
		dd <- density(xx[xx < Th], na.rm=TRUE, ...)
		# bg <- dd$x[which.max(dd$y)]
		bg <- switch(bgMethod,
			none = 0,
			density = dd$x[which.max(dd$y)],
			mean = mean(xx),
			median = median(xx))
		bgNum <- 2 * length(which(xx <= bg))
		# estimate the foreground level by removing the background probes
		fg <- switch(fgMethod,
			density = (sum(xx) - bgNum * bg) / (length(xx) - bgNum),
			mean = mean(xx),
			median = median(xx))
		twoPoint <- c(bg, fg)
		return(twoPoint)
	})

	if (is.null(targetArray)) {
		target2P <- apply(twoPoint, 1, mean)
	} else {
		target2P <- twoPoint[, targetArray]
	}
	
	## normalization based on twoPoint (background and foreground levels)
	normalized <- expr
	nArray <- ncol(expr)
	for (i in 1:nArray) {
		if (scaling) {
			ff <- abs((target2P[2] - target2P[1]) / (twoPoint[2,i] - twoPoint[1,i]))
		} else {
			ff <- 1
		}
		normalized[,i] <- (expr[,i] - twoPoint[1,i]) * ff + target2P[1]
	}

	## if the targetArray is an external vector, it will be removed from the normalized data.
	if (externalTarget) normalized <- normalized[,-1]
	## transformed as original scale in not log2transformed
	if (log2Trans) {
		if (min(normalized) <= 0) {
			normalized <- lumiB(normalized, method='forcePositive')
		}
		normalized <- log2(normalized)	
	} 
	
	if (is(x.lumi, 'ExpressionSet')) {
		exprs(x.lumi) <- normalized
	} else {
		x.lumi <- normalized
	}
	
	return(x.lumi)
}
