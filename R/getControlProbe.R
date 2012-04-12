`getControlProbe` <-
function(controlData, type=NULL) 
{
	if (is(controlData, 'LumiBatch')) {
		controlData <- controlData@controlData
	} 
	allControlType <- controlData$controlType
	uniControlType <- getControlType(controlData)
	allProbeID <- controlData$ProbeID
	selProbeID <- allProbeID
	if (!is.null(type)) {
		type <- toupper(type)
		allControlType <- toupper(allControlType)
		if (type %in% uniControlType) {
			selProbeID <- allProbeID[allControlType == type[1]]
		} 
	} 
	return(selProbeID)
}

