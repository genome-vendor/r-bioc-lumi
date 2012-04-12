`monoSmu` <-
function(x, y, newX=NULL, 
	nSupport=min(200, length(x)),   	# downsampled data points
	nKnots=6, rotate=FALSE, ifPlot=FALSE, xlab='x', ylab='y', ...
	) {
	
	x.old <- x
	y.old <- y
	
	if (rotate) {
		y <- y.old - x.old		## Equivalent to M
		x <- y.old + x.old 		## Equivalent to 2*A
	}
	
	# the order of (x,y) pairs does not matter
	# order of newX does matter
	ord <- order(x)
	x.ord <- x[ord]; y.ord <- y[ord]
	
	xy <- data.frame(x=x.ord, y=y.ord)
	result.sm <- supsmu(x=xy$x, y=xy$y, ...)
	
	if (rotate) {
		## Rotate back to the original direction
		x.temp <- (result.sm$x - result.sm$y)/2
		y.temp <- (result.sm$x + result.sm$y)/2
		xy1 <- data.frame(x=x.temp, y=y.temp)
	} else {
	    xy1 <- data.frame(x=result.sm$x, y=result.sm$y)
	}
	
	## down-sampling
	rangeX <- range(xy1$x)
	selInd <- round(seq(from=1, to = nrow(xy1), length =nSupport))
	## make sure there are enough samples at high end
	selX <- seq(rangeX[1] + diff(rangeX)/2, rangeX[2], length=10)
	selIndHigh <- sapply(selX, function(x) which.min(abs(xy1$x - x)))
	selInd <- sort(c(selInd, selIndHigh))
	xy1 <- xy1[selInd, ]
	
	if (is.null(newX))  newX <- x.old 
	
	newY <- monoSpline(x=xy1$x, y=xy1$y, newX=newX, nKnots=nKnots, ifPlot=FALSE)

	if (ifPlot) {
		plot(x.old, y.old, pch='.', xlab=xlab, ylab=ylab)
		points(xy1$x, xy1$y, col=2)
		ii <- order(newX)
		lines(newX[ii], newY[ii], col=3, type='l')
	}
	
	return(newY) 
}

