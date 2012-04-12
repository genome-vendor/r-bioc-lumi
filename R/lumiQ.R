`lumiQ` <-
function(x.lumi, logMode=TRUE, detectionTh=0.01) {
	if (!is(x.lumi, 'LumiBatch')) stop('The object should be class "LumiBatch"!')
	history.submitted <- as.character(Sys.time())

	expr <- exprs(x.lumi)
	if (any(is.na(expr))) {
		naInd <- apply(expr, 1, function(x) any(is.na(x)))
		expr <- expr[!naInd,]
	}

	if (logMode && (max(expr, na.rm=TRUE) > 50)) {
		# remove the negative values
		if (min(expr) < 0) {
			rMin <- rowMin(expr)
			expr <- expr[rMin > 0, ]
		}
		expr <- log2(expr)
	} 
	sampleName <- colnames(expr)

	## mean, variance, 'cv', 'AP', 'density', 'correlation', 'sample relation'
	mm <- colMeans(expr)
	std <- apply(expr, 2, sd)
	
	## AP calls
	if (!is.null(x.lumi@QC$sampleSummary)) {
		detectionRate <- x.lumi@QC$sampleSummary[3,]
		detectionName <- rownames(x.lumi@QC$sampleSummary)[3]
	} else {
		detectionName <- paste('detection rate(', detectionTh, ')', sep='')
		if (!is.null(detection(x.lumi))) {
			detectionRate <- detectionCall(x.lumi, Th=detectionTh, type='sample') / nrow(x.lumi)
		} else {
			detectionRate <- rep(NA, ncol(x.lumi))
		}
	}

	## detect outlier
	center <- rowMeans(expr)
	profile <- cbind(center, expr)
	colnames(profile) <- c('Center', colnames(expr))
	distCenter <- as.matrix(dist(t(profile), method="euclidean"))

	sampleSummary <- rbind(mm, std, detectionRate, distCenter[2:nrow(distCenter),1])
	rownames(sampleSummary) <- c('mean', 'standard deviation', detectionName, 'distance to sample mean')
	sampleSummary <- signif(sampleSummary, 4)

	## record history
	history.finished <- as.character(Sys.time())
	history.command <- capture.output(print(match.call(lumiQ)))

	if (is.null(x.lumi@history$lumiVersion)) x.lumi@history$lumiVersion <- rep(NA, nrow(x.lumi@history))
	lumiVersion <- packageDescription('lumi')$Version
	x.lumi@history <- rbind(x.lumi@history, data.frame(submitted=history.submitted, 
			finished=history.finished, command=history.command, lumiVersion=lumiVersion))

	## create a QC slot
	QC <- x.lumi@QC
	QC$sampleSummary <- sampleSummary
	QC$history <- x.lumi@history
	
	x.lumi@QC <- QC
	return(x.lumi)
}

