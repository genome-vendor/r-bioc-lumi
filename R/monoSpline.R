`monoSpline` <-
function(x, y, newX=NULL, nKnots=6, ifPlot=FALSE) {
	# behavior: if extropolation, outside values are constant to the
	#  nearest inside one
	if(!require(mgcv)) stop('Package "mgcv" should be installed!')
	
	nKnots <- round(nKnots)
	if (is.null(newX)) newX <- x
	
	## make sure there is no other "s" functions before mgcv::s in the search path
	rmPackage <- find('s')
	if (rmPackage[1] != "package:mgcv") {
		detach("package:mgcv")
		library(mgcv)
	}
	
	# Show regular spline fit (and save fitted object)
	f.ug <- gam(y ~ s(x, k=nKnots, bs="cr"))
	
	# Create Design matrix, constraints etc. for monotonic spline....
	dat <- data.frame(x=x, y=y)
	sm <- smoothCon(s(x, k=nKnots, bs="cr"), dat, knots=NULL)[[1]]
	if (length(sm$xp) < 6) warning('Few than 6 nKnots were specified!\n Possible inaccurate fitting!\n')
	F <- mono.con(sm$xp);   # get constraints
	G <- list(X=sm$X, C=matrix(0,0,0), sp=f.ug$sp, p=sm$xp, y=y, w=y*0+1, 
			Ain=F$A, bin=F$b, S=sm$S, off=0)
	
	p <- pcls(G);  # fit spline(using s.p. from unconstrained fit)
	fv <- Predict.matrix(sm, data.frame(x=newX)) %*% p
	fv <- as.vector(fv)
	
	if (ifPlot) {
	    plot(x,y); lines(newX, fv, col=2)
	}
	
	return(fv)
}

