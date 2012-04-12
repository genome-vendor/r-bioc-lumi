`addControlData2lumi` <-
function(controlData, x.lumi) 
{
	if (missing(x.lumi) || missing(controlData)) stop('Both controlData and x.lumi are required!')
	if (is.character(controlData)) {
		controlFile <- controlData
		allControlInfo <- lumiR.batch(controlFile, lib=NULL, checkDupId=FALSE)
		controlData <- as.data.frame(exprs(allControlInfo))
		controlType <- as.character(pData(featureData(allControlInfo))$TargetID)
		ProbeID <- as.character(pData(featureData(allControlInfo))$ProbeID)
		controlData <- data.frame(controlType=controlType, ProbeID=ProbeID, controlData)
	}
	if (is.matrix(controlData)) controlData <- as.data.frame(controlData)
	if (is(controlData, 'data.frame')) {
		## match the column names of controlData and LumiBatch object
		sampleID <- as.character(pData(phenoData(x.lumi))$sampleID)
		if (is.null(sampleID)) sampleID <- sampleNames(x.lumi)
		controlSampleID <- names(controlData)
		if ('TargetID' %in% controlSampleID) {
			controlSampleID[controlSampleID == 'TargetID'] <- 'controlType'
		}
		if (all(sampleID %in% controlSampleID)) {
			x.lumi@controlData <- controlData[, c('controlType', 'ProbeID', sampleID)]
		} else {
			sampleIDInfo <- strsplit(sampleID, split="_")
			newID <- NULL
			temp <- lapply(sampleIDInfo, function(x) {
				newID <<- c(newID, paste(x[1:2], collapse="_"))
			})
			if (all(newID %in% controlSampleID)) {
				x.lumi@controlData <- controlData[, c('controlType', 'ProbeID', newID)]
			} else {
				stop('SampleID does not match up between controlData and x.lumi!')				
			}
		}
		names(x.lumi@controlData) <- c('controlType', 'ProbeID', sampleNames(x.lumi))		
	} else {
		stop('Input data type is not supported!')
	}
	return(x.lumi)
}

