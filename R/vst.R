`vst` <-
function(u, std, nSupport=min(length(u), 500), backgroundStd=NULL, fitMethod=c('linear', 'quadratic'), lowCutoff=1/3, ifPlot=FALSE) {
	# u is the mean of probe beads
	# std is the standard deviation of the probe beads
	# 
	fitMethod <- match.arg(fitMethod)
	## Estimate the background variance c3
	c3 <- ifelse (is.null(backgroundStd), 0, backgroundStd)
	
	ord <- order(u); u.bak <- u
	u <- u[ord]; std <- std[ord]
	
	if (any(std < 0)) {
		stop('Negative expression standard deviation is not allowed!')
	}
	

	if (fitMethod == 'quadratic') {
		downSampledU <- u; offset=0
		dd <- data.frame(y=std, x2=u^2, x1=u)
		lm2 <- lm(y ~ x2 + x1, dd)
		smoothStd <- predict(lm2, data.frame(x2=downSampledU^2, x1=downSampledU))
	} else {
		## downsampling to speed up 
		if (min(u) < 1) {
			offset <- 1 - min(u)
		} else {
			offset <- 0
		}
		offset <- 1 - min(u)
		downSampledU <- 2^seq(from=log2(min(u + offset)), to=log2(max(u + offset)), length=nSupport) - offset

		minU <- log2(max(100 + offset, min(u)))
		maxU <- log2(max(u))
		# uCutoff <- 2^((maxU + minU)/2)
		uCutoffLow <- 2^(minU + (maxU - minU) * lowCutoff)
		uCutoffHigh <- 2^(minU + (maxU - minU) * 4/5)
		selInd <- (u > uCutoffLow & u < uCutoffHigh)
		selLowInd <- (u < uCutoffLow)
		if (c3 != 0) {
			selInd <- selInd & (std^2 > c3)
			dd <- data.frame(y=sqrt(std[selInd]^2 - c3), x1=u[selInd])
 			# if (nrow(dd) > 5000) dd <- dd[sample(1:nrow(dd), 5000),]
			lmm <- lm(y ~ x1, dd)
			c1 <- lmm$coef[2]
			c2 <- lmm$coef[1]
		} else {
			iterNum <- 0
			c3.i <- 0
			while(iterNum <= 20) {
				selInd.i <- selInd & (std^2 > c3.i)
				dd <- data.frame(y=sqrt(std[selInd.i]^2 - c3.i), x1=u[selInd.i])
				# if (nrow(dd) > 5000) dd <- dd[sample(1:nrow(dd), 5000),]
				lm.i <- lm(y ~ x1, dd)
				c1.i <- lm.i$coef[2]
				c2.i <- lm.i$coef[1]
				y <- std[selLowInd]
				x <- u[selLowInd]
				cc <- y^2 - (c1.i * x + c2.i)^2
				c3.i.new <- mean(cc, trim=0.05)
				if (c3.i.new < 0) {
					break
				} else {
					if (abs(c3.i.new - c3.i) < 1e-5) break
					c3.i <- c3.i.new
				}
				iterNum <- iterNum + 1
			}
			c1 <- c1.i; c2 <- c2.i; c3 <- c3.i
			if (c3 < 0) c3 <- 0
		}
		smoothStd <- ((c1 * downSampledU + c2)^2 + c3)^(1/2)
	}

	if (ifPlot) {
		par(mfrow=c(1,2))
		x <- u[u > 0]
		y <- std[u > 0]
		len <- length(x)
		ind <- sample(1:len, min(5000, len))
		plot(x[ind], y[ind], pch='.', log='xy', xlab="Mean", ylab="Standard Deviation", main='(A) Relations of probe Mean and SD')
		lines(downSampledU, smoothStd, col=3, lwd=1.5)
	}
	
	## calculate the integration (h function is the integral)
	if (fitMethod == 'linear') {
		if (c3 == 0) {
			## Transform function h(x) = g * log(a + b * x)
			g <- 1/c1
			a <- c2
			b <- c1
			tmp <- a + b * u.bak
			if (any(tmp < 0)) {
				transformedU <- log(u.bak)
				g <- 1; a <- 0; b <- 1
			} else {
				transformedU <- g * log(a + b * u.bak)
			}
			if (ifPlot) hy <- g * log(a + b * downSampledU) 
			transFun <- 'log'
		} else {
			## Transform function h(x) = g * asinh(a + b * x)
			g <- 1/c1
			a <- c2/sqrt(c3)
			b <- c1/sqrt(c3)
			transformedU <- g * asinh(a + b * u.bak)
			if (ifPlot) hy <- g * asinh(a + b * downSampledU) 
			transFun <- 'asinh'
		}
		transform.parameter <- c(a, b, g, 0)
		names(transform.parameter) <- c('a', 'b', 'g', 'Intercept')
	} else {
		dd <- diff(downSampledU)	## the interval between samples
		dd <- c(dd[1], dd)
		hy <- cumsum(1/smoothStd * dd)		
		# get a smoothed version and interpolation of original data order
	    transformedU <- monoSpline(x=downSampledU, y=hy, newX=u.bak, nKnots=20, ifPlot=FALSE)
		transFun <- 'quadratic'
		transform.parameter <- NULL
	}	

	## rescale to the similar range with log2
	if (fitMethod == 'quadratic') {
		minU <- log2(max(100 + offset, min(u)))
		maxU <- log2(max(u))
		uCutoffLow <- 2^(minU + (maxU - minU) * lowCutoff)
		uCutoffHigh <- 2^(minU + (maxU - minU) * 4/5)		
	}
	cutInd <- which.min(abs(u.bak - uCutoffLow))
	maxInd <- which.max(u.bak)
	y <- c(u.bak[cutInd], u.bak[maxInd])
	x <- c(transformedU[cutInd], transformedU[maxInd])
	m <- lm(log2(y) ~ x)
	if (fitMethod == 'linear') {
		transform.parameter <- c(a, b, g * m$coef[2], m$coef[1])
		names(transform.parameter) <- c('a', 'b', 'g', 'Intercept')
		## The transform parameter is in the transFun below
		## transFun <- g * asinh(a + b * x) * m$coef[2] + m$coef[1]
	} else {
		transform.parameter <- NULL
		transFun <- NULL
	}
	transformedU <- predict(m, data.frame(x=transformedU))
	attr(transformedU, 'parameter') <- transform.parameter
	attr(transformedU, 'transformFun') <- transFun

    if (ifPlot) {
		x <- log2(u.bak[u.bak > 0])
		y <- transformedU[u.bak > 0]
        ii <- order(x)
        plot(x[ii], y[ii], lwd=1.5, type='l', xlab="Log2 transformed value", ylab="VST transformed value", main='(B) Log2 vs. VST')
		abline(a=0, b=1, col=3, lwd=1.5, lty=2)
		par(mfrow=c(1,1))
    }
	
	return(transformedU)
}

