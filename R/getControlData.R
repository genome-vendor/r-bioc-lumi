`getControlData` <-
function(x, type=c('data.frame', 'LumiBatch')) 
{
	type <- match.arg(type)
	if (is.character(x)) {
		allControlInfo <- lumiR.batch(x, lib=NULL, checkDupId=FALSE, convertNuID=FALSE)
		if (type == 'LumiBatch') {
			return(allControlInfo)
		}
		x <- allControlInfo
	} 

	if (is(x, 'LumiBatch')) {
		if (type == 'LumiBatch') {
			return(x)
		} else {
			if (nrow(x@controlData) == 0) {
				controlData <- as.data.frame(exprs(x))
				controlType <- pData(featureData(allControlInfo))$TargetID
				if (length(which(toupper(controlType) == 'NEGATIVE')) > 10) {
					ProbeID <- pData(featureData(allControlInfo))$ProbeID
					controlNames <- names(controlData)
					controlData <- data.frame(controlType=as.character(controlType), ProbeID=as.character(ProbeID), controlData)
					names(controlData) <- c('controlType', 'ProbeID', controlNames)
				} else {
					controlData <- data.frame()
				}
			} else {
				controlData <- x@controlData
			}
		}		
	} else {
		stop('Input data should be a control data file or a LumiBatch object!')
	}
	return(controlData)
}

