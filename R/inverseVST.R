`inverseVST` <-
function(x, fun=c('asinh', 'log'), parameter=NULL) {
	if (missing(x)) {
		stop('Please provide the data!')
	}
	if (!is.null(attr(x, 'vstParameter'))) {
		parameter <- attr(x, 'vstParameter')
		fun <- attr(x, 'transformFun')
	} else {
		if (is.null(parameter)) stop('Please provide the transform parameters!')
		fun <- match.arg(fun)
		## parameter is in the format of: "a, b, g, Intercept"
		a <- parameter[1]
		b <- parameter[2]
		g <- parameter[3]
		intercept <- parameter[4]
	}
	if (is(x, 'ExpressionSet')) {
		inv <- x
		attr(inv, 'vstParameter') <- NULL
		attr(inv, 'transformFun') <- NULL
		if (length(fun) == 1) {
			a <- parameter[1]
			b <- parameter[2]
			g <- parameter[3]
			intercept <- parameter[4]
			if (fun == 'asinh') {
				inv.x <- (sinh((exprs(x) - intercept)/g) - a)/b
			} else if (fun == 'log') {
				inv.x <- (exp((exprs(x) - intercept)/g) - a)/b
			}
			rownames(inv.x) <- rownames(exprs(x))
			colnames(inv.x) <- colnames(exprs(x))
			exprs(inv) <- inv.x
		} else {
			if (length(fun) != ncol(inv)) stop('The length of "fun" does not match the data!')
			for (i in 1:length(fun)) {
				fun.i <- fun[i]
				parameter.i <- parameter[i,]
				a <- parameter.i[1]
				b <- parameter.i[2]
				g <- parameter.i[3]
				intercept <- parameter.i[4]
				if (fun.i == 'asinh') {
					exprs(inv)[,i] <- (sinh((exprs(x)[,i] - intercept)/g) - a)/b
				} else if (fun.i == 'log') {
					exprs(inv)[,i] <- (exp((exprs(x)[,i] - intercept)/g) - a)/b
				}
			}
		}
	} else {
		if (fun == 'asinh') {
			## Transform function h(x) = g * asinh(a + b * x)
			## transFun <- g * asinh(a + b * x) * m$coef[2] + m$coef[1]
			inv <- (sinh((x - intercept)/g) - a)/b
		} else if (fun == 'log') {
			## Transform function h(x) = g * log(a + b * x)
			## transFun <- g * log(a + b * x) * m$coef[2] + m$coef[1]
			inv <- (exp((x - intercept)/g) - a)/b
		}
	}
	
	return(inv)
}

