`plotControlData` <-
function(controlData, type=NULL, slideIndex=NULL, logMode=FALSE, new=TRUE, ...) 
{
	if (is(controlData, 'LumiBatch')) {
		sampleID <- pData(phenoData(controlData))$sampleID
		controlData <- controlData@controlData
		if (nrow(controlData) == 0) stop('Slot controlData is empty!')
	} else {
		sampleID <- NULL
	}
	allControlType <- controlData$controlType
	uniControlType <- getControlType(controlData)
	controlData <- controlData[, -c(1,2)]
	if (!is.null(type)) {
		type <- toupper(type)
		allControlType <- toupper(allControlType)
		if (!all(type %in% uniControlType)) {
			warning('"type" does not match "controlType". Please use getControlType function to view the available type!')
			type <- ''
		} 
	} else {
		type <- ''
	}
	if (is.null(slideIndex)) {
		if (is.null(sampleID)) sampleID <- colnames(controlData)
		sampleIDInfo <- strsplit(sampleID, split="_")
		chipID <- sapply(sampleIDInfo, function(x) x[1])
		ord <- order(chipID)
		chipID <- chipID[ord]
		controlData <- controlData[, ord]
		col <- as.numeric(as.factor(chipID))
		chipNum <- length(unique(chipID))
	} else {
		ord <- order(slideIndex)
		controlData <- controlData[, ord]
		slideIndex <- slideIndex[ord]
		col <- as.numeric(as.factor(slideIndex))
		chipNum <- length(unique(slideIndex))
	}
	if (logMode) {
		if (max(controlData) > 50) {
			controlData <- controlData - min(controlData) + 1
			controlData <- log2(controlData)
		}
		ylab <- 'Expression Amplitude (log2)'
	} else {
		if (max(controlData) < 50) controlData <- 2^(controlData)
		ylab <- 'Expression Amplitude'
	}
	
	if (length(type) > 1) oldPar <- par(mfrow=c(ceiling(length(type)/2), 2))
	for (type.i in type) {
		if (type.i != '') {
			selControlData <- controlData[allControlType == type.i, , drop=FALSE]
		} else {
			selControlData <- controlData
		}
		
		if (nrow(selControlData) > 1) {
			mm <- apply(selControlData, 2, mean)
			std <- apply(selControlData, 2, sd)
		} else {
			mm <- selControlData
			std <- rep(0, length(mm))
		}

		labels <- colnames(selControlData)
		if (is.null(labels)) labels <- as.character(1:ncol(selControlData))
		## set the margin of the plot
		mar <- c(max(nchar(labels))/2 + 4.5, 5, 5, 3)
		old.mar <- par('mar')
		old.xaxt <- par('xaxt')
		par(xaxt='n')
		par(mar=mar)
		
		ind <- ncol(selControlData)
		if (new) {
			range <- c(min(mm - std * 1.5), max(mm + std * 1.1))
			plot(1:ind, mm, pch = 19, cex=1.8, type='o', col=col, xlab='', ylab='Expression Amplitude', ylim=range, ...)
			title(type.i)
		} else {
			points(1:ind, mm, pch = 19, cex=1.8, type='o', col=col, ...)
		}
		if (nrow(selControlData) > 1) {
			arrows(1:ind, mm - std, 1:ind, mm + std, code = 3, col=col, angle = 90, length = .1)
		}
		# plot the vertical lines
		if (chipNum > 1) {
			abline(v=0.5 + which(diff(col) != 0), lty=2)
		}

		par(xaxt='s')
		axis(1, at=1:ncol(selControlData), labels=labels, tick=TRUE, las=2)
		par(mar=old.mar)
		par(xaxt=old.xaxt)
	}
	if (length(type) > 1) par(oldPar)
	
	return(invisible(TRUE))
}

