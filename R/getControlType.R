`getControlType` <-
function(controlData) 
{
	if (is(controlData, 'LumiBatch')) {
		controlData <- controlData@controlData
	} 
	if (is(controlData, 'data.frame')) {
		if (nrow(controlData) == 0) stop('controlData is empty!')
		return(unique(as.character(controlData$controlType)))
	} else {
		return(NA)
	}
}

