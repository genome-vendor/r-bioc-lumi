`detectOutlier` <-
function(x, metric="euclidean", standardize=TRUE, Th=2, ifPlot=FALSE) {
	
	if (is(x, 'ExpressionSet')) {
	    profile <- exprs(x)		
	} else if (is.numeric(x)) {
		profile <- as.matrix(x)
	} else {
		stop('The object should be a matrix or class "ExpressionSet" inherited!')
	}
	
	center <- rowMeans(profile)
	profile1 <- cbind(center, profile)
	colnames(profile1) <- c('Center', colnames(profile))


	if (metric == 'cor') {
		d <- cor(profile1)
		mad.d <- median(1 - d[2:nrow(d),1])
	} else {
		if (standardize) {
			profile1 <- scale(profile1)
		}
		d <- as.matrix(dist(t(profile1), method=metric))
		mad.d <- median(d[2:nrow(d),1])
	}

	Th <- Th * mad.d
	outlier <- (d[2:nrow(d),1] >= Th)
	attr(outlier, 'sampleDistance') <- d
	attr(outlier, 'threshold') <- Th
	main <- paste('Outlier detection based on sample distance to "Center"')

	if (ifPlot) {
		hc = hclust(as.dist(d), 'ave')
		plot(hc, xlab='Sample', ylab='Distance', main=main)
		abline(h=Th, col=2, lty=2)
		return(invisible(TRUE))	
	} else {
		return(outlier)
	}
}

