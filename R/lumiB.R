lumiB <- function(x.lumi, method = c('none', 'bgAdjust', 'forcePositive', 'bgAdjust.affy'), ...) 
{

	if (is(x.lumi, 'ExpressionSet')) {
	    # x.lumi is a LumiBatch object
	    x.matrix <- exprs(x.lumi)		
	} else if (is.numeric(x.lumi)) {
		x.matrix <- as.matrix(x.lumi)
	} else {
		stop('The object should be a matrix or class "ExpressionSet" inherited!')
	}
	
	if (!(is.function(method)) && !(method %in% c('bgAdjust', 'none', 'forcePositive', 'bgAdjust.affy'))) {
		print('This method is not supported!')
		return(x.lumi)
	} else if (method == 'none') {
		return(x.lumi)
	} else if (method == 'forcePositive') {
		if (min(x.matrix) > 0)  return(x.lumi)
	}
	
	history.submitted <- as.character(Sys.time())
	if (method == 'bgAdjust') {
		return(bgAdjust(x.lumi, ...))
	} else if (method == 'bgAdjust.affy') {
		x.matrix <- apply(x.matrix, 2, bg.adjust, ...) 
	} else if (method == 'forcePositive') {
		offset <- apply(x.matrix, 2, min)
		offset[offset <= 0] <- offset[offset <= 0] - 1.01
		offset[offset > 0] <- 0
		offset <- rep(1, nrow(x.matrix)) %*% t(offset)
		x.matrix <- x.matrix - offset
	} else if (is.function(method)) {
		x.lumi <- method(x.lumi, ...)
	} else {
		print('The method is not supported!')
		return(x.lumi)
	}

	if (is(x.lumi, 'ExpressionSet')) {
		if (!is.function(method)) exprs(x.lumi) <- x.matrix
		if (is(x.lumi, 'LumiBatch')) {
			# history tracking
			history.finished <- as.character(Sys.time())
			history.command <- capture.output(print(match.call(lumiB)))
        	
			if (is.null(x.lumi@history$lumiVersion)) x.lumi@history$lumiVersion <- rep(NA, nrow(x.lumi@history))
			lumiVersion <- packageDescription('lumi')$Version
			x.lumi@history<- rbind(x.lumi@history, data.frame(submitted=history.submitted, 
					finished=history.finished, command=history.command, lumiVersion=lumiVersion))
		}
	} else {
		x.lumi <- x.matrix
	}
	return(x.lumi)
}