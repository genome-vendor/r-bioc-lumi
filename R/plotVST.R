`plotVST` <-
function(x, transFun = NULL, plotRange = NULL, addLegend = TRUE, ...) {
	if (missing(x)) {
		stop('Please provide all input parameters!')
	} else {
		parameter <- NULL
		if (is(x, 'LumiBatch')) {
			parameter <- attr(x, 'vstParameter')
			transFun <- attr(x, 'transformFun')
			if (is.null(plotRange)) {
				plotRange <- range(exprs(x))
				if (max(plotRange) < 100) plotRange <- 2^plotRange
			}
		}
		if (is.null(plotRange)) plotRange <- c(10, 10000)
		if (is(x, 'data.frame')) parameter <- as.matrix(x)
		if (length(x) == 4) {
			parameter <- matrix(x, nrow=1)
			colnames(parameter) <- names(x)
		}
		## parameter is in the format of: "a, b, g, Intercept"
		if (!all(colnames(parameter) %in% c('a', 'b', 'g', 'Intercept')))
			stop('Incorrect format of input x!')
		if (is.null(transFun)) transFun <- rep('asinh', nrow(parameter))
		if (length(transFun) != nrow(parameter)) stop('Length of "transFun" does not match "x"!')
	}
	
	xx <- NULL
	trans <- NULL
	for (i in 1:nrow(parameter)) {
		a <- parameter[i, 'a']
		b <- parameter[i, 'b']
		g <- parameter[i, 'g']
		intercept <- parameter[i, 'Intercept']
		fun <- transFun[i]
		xx.i <- seq(plotRange[1], plotRange[2], length=500)
		if (fun == 'asinh') {
			## Transform function h(x) = g * asinh(a + b * x)
			## transFun <- g * asinh(a + b * x) * m$coef[2] + m$coef[1]
			trans.i <- g * asinh(a + b * xx.i) + intercept
		} else if (fun == 'log') {
			## Transform function h(x) = g * log(a + b * x)
			## transFun <- g * log(a + b * x) * m$coef[2] + m$coef[1]
			trans.i <- g * log(a + b * xx.i) + intercept
		}
		xx <- cbind(xx, xx.i)
		trans <- cbind(trans,trans.i)
	}
	matplot(xx, trans, type='l', lty=1:nrow(parameter), col=1:nrow(parameter), xlab='un-transformed value', ylab='transformed value', ...)
	if (addLegend) {
		r.trans <- range(trans)
		legend(x=max(xx)*2/3, y=min(r.trans) + diff(r.trans)/3, legend=rownames(parameter), lty=1:nrow(parameter), col=1:nrow(parameter))
	}
	return(invisible(list(untransformed=xx, transformed=trans)))
}

