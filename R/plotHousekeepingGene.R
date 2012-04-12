`plotHousekeepingGene` <-
function(controlData, lib=NULL, slideIndex=NULL, addLegend=TRUE, logMode=TRUE, ...) 
{
	if (is(controlData, 'LumiBatch')) {
		sampleID <- pData(phenoData(controlData))$sampleID
		controlData <- controlData@controlData
		if (nrow(controlData) == 0) stop('Slot controlData is empty!')
	} else {
		sampleID <- NULL
	}
	allControlType <- controlData$controlType
	allControlProbe <- controlData$ProbeID
	uniControlType <- getControlType(controlData)
	controlData <- controlData[, -c(1,2)]

	if (is.null(slideIndex)) {
		if (is.null(sampleID)) sampleID <- colnames(controlData)
		sampleIDInfo <- strsplit(sampleID, split="_")
		chipID <- sapply(sampleIDInfo, function(x) x[1])
		ord <- order(chipID)
		chipID <- chipID[ord]
		controlData <- controlData[, ord]
		chip <- as.numeric(as.factor(chipID))
		chipNum <- length(unique(chipID))
	} else {
		ord <- order(slideIndex)
		controlData <- controlData[, ord]
		slideIndex <- slideIndex[ord]
		chip <- as.numeric(as.factor(slideIndex))
		chipNum <- length(unique(slideIndex))
	}
	
	# get HOUSEKEEPING
	ind <- grep('house', allControlType, ignore.case=TRUE)
	if (length(ind) == 0) stop('No housekeeping gene information found!')
	selControlData <- controlData[ind, ]
	hkgene <- allControlProbe[ind]
	if (!is.null(lib)) {
		require(annotate)
		require(lib, character.only=TRUE)
		hkSymbol <- getSYMBOL(hkgene, lib)
		hkSymbol[is.na(hkSymbol)] <- hkgene[is.na(hkSymbol)]
		hkgene <- hkSymbol
	}
	
	if (logMode) {
		if (max(selControlData) > 50) {
			selControlData <- selControlData - min(selControlData) + 1
			selControlData <- log2(selControlData)
		}
		ylab <- 'Expression Amplitude (log2)'
	} else {
		if (max(selControlData) < 50) selControlData <- 2^(selControlData)
		ylab <- 'Expression Amplitude'
	}
	labels <- colnames(selControlData)
	if (is.null(labels)) labels <- as.character(1:ncol(selControlData))
	## set the margin of the plot
	mar <- c(max(nchar(labels))/2 + 4.5, 5, 5, 3)
	old.mar <- par('mar')
	old.xaxt <- par('xaxt')
	par(xaxt='n')
	par(mar=mar)
	
	col <- lty <- 1:length(hkgene)
	if (addLegend) {
		xlim <- c(1, ncol(selControlData) + 1)
	} else {
		xlim <- c(1, ncol(selControlData))
	}
	matplot(t(selControlData), pch = 19, cex=1.2, type='o', col=col, lty=lty, xlim=xlim,
	 	xlab='', ylab=ylab)
	title('Expression profile of housekeeping genes')
	# plot the vertical lines
	if (chipNum > 1) abline(v = 0.5 + which(diff(chip) != 0), lty=2)
	if (addLegend) legend('topright', legend=hkgene, lty=lty, col=col, inset = .02)

	par(xaxt='s')
	axis(1, at=1:ncol(selControlData), labels=labels, tick=TRUE, las=2)
	par(mar=old.mar)
	par(xaxt=old.xaxt)
	
	return(invisible(TRUE))
}

