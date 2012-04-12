##=================================================
## Define LumiBatch object:

setClass('LumiBatch', 
	representation(history='data.frame', controlData='data.frame', QC='list'), 
	prototype=list(history=data.frame(
		submitted   = I(vector()),
		finished    = I(vector()),
		command     = I(vector()),
		lumiVersion = I(vector())
	), controlData = data.frame(), QC = list()),
	contains='ExpressionSet')


setMethod('initialize', 'LumiBatch', function(.Object, 
	exprs = new('matrix'),
	se.exprs = new('matrix'),		# standard deviation of the bead measurements of each probe
	...) 
{
	callNextMethod(.Object, exprs=exprs, se.exprs=se.exprs, ...)
})


setValidity("LumiBatch", function(object) 
{
    msg <- Biobase:::validMsg(NULL, Biobase:::isValidVersion(object, "ExpressionSet"))
    msg <- Biobase:::validMsg(msg, assayDataValidMembers(assayData(object), c("exprs")))
    msg <- Biobase:::validMsg(msg, assayDataValidMembers(assayData(object), c("se.exprs")))
    if (is.null(msg)) TRUE else msg
})


##=================================================
## methods
if (is.null(getGeneric("getHistory"))) setGeneric("getHistory", function(object) standardGeneric("getHistory"))
if (is.null(getGeneric("beadNum"))) setGeneric("beadNum", function(object) standardGeneric("beadNum"))
if (is.null(getGeneric("beadNum<-"))) setGeneric("beadNum<-", function(object, value) standardGeneric("beadNum<-"))

if (is.null(getGeneric("detection"))) setGeneric("detection", function(object) standardGeneric("detection"))
if (is.null(getGeneric("detection<-"))) setGeneric("detection<-", function(object, value) standardGeneric("detection<-"))

if (is.null(getGeneric("summary"))) setGeneric("summary", function(object, ...) standardGeneric("summary"))
if (is.null(getGeneric("show"))) setGeneric("show", function(object) standardGeneric("show"))
if (is.null(getGeneric("combine"))) setGeneric("combine", function(x, y, ...) standardGeneric("combine"))
if (is.null(getGeneric("density"))) setGeneric("density", function(x, ...) standardGeneric("density"))

setMethod("se.exprs", signature(object="LumiBatch"),
          function(object) assayDataElement(object,"se.exprs"))

setReplaceMethod("se.exprs", signature(object="LumiBatch",value="matrix"),
                 function(object, value) assayDataElementReplace(object, "se.exprs", value))

setMethod("beadNum", signature(object="LumiBatch"),
          function(object) assayDataElement(object,"beadNum"))

setReplaceMethod("beadNum", signature(object="LumiBatch",value="matrix"),
                 function(object, value) assayDataElementReplace(object, "beadNum", value))

setMethod("detection", signature(object="LumiBatch"),
          function(object) assayDataElement(object,"detection"))

setReplaceMethod("detection", signature(object="LumiBatch",value="matrix"),
                 function(object, value) assayDataElementReplace(object, "detection", value))

setMethod("getHistory",signature(object="LumiBatch"), function(object) object@history)


setMethod("summary",signature(object="LumiBatch"), function(object, type=c('data', 'QC')) 
{
	type <- match.arg(type)
	if (type == 'data') {
		show(object)		
	} else {
		if (is.null(object@QC$sampleSummary)) {
			cat('Run Quality Control estimation ...\n')
			object <- lumiQ(object)
		}
		if (!is.null(object@QC$sampleSummary)) {
			QC <- object@QC
			dimen <- dim(exprs(object))
			cat(paste("Data dimension: ", paste(dimen[1], 'genes', 'x', dimen[2], 'samples', collapse=" "), '\n'))
			sampleSummary <- QC$sampleSummary
			cat(paste("\nSummary of Samples:\n", sep=''))
			print(sampleSummary, quote=FALSE)

			cat('\nMajor Operation History:\n')
			print(QC$history, quote=FALSE) 
		}
	}
})


setMethod("show",signature(object="LumiBatch"), function(object) 
{
	cat('Summary of BeadStudio output:\n\t')
	cat(notes(object)[[1]], sep='\n\t')
	cat('Major Operation History:\n')
	print(getHistory(object)) 
	cat('\nObject Information:\n')
	callNextMethod()
	if (!is.null(object@controlData)) {
		cat('Control Data: Available\n')
	} else {
		cat('Control Data: N/A\n')
	}
	cat("QC information: Please run summary(x, 'QC') for details!\n")
})


setMethod("[", "LumiBatch", function(x, i, j, ..., drop = FALSE) 
{
	if (missing(drop)) drop <- FALSE
   	history.submitted <- as.character(Sys.time())
	ddim <- dim(x)
	
	sampleName <- sampleNames(x)
	## do default processing of 'ExpressionSet'
	x <- callNextMethod()

	if (!missing(i) && !missing(j)) {
		history.command <- paste('Subsetting', ddim[1], 'features and', ddim[2], 'samples.')		
	} else if (!missing(i)) {
		history.command <- paste('Subsetting', ddim[1], 'features.')
	} else if (!missing(j)) {
		history.command <- paste('Subsetting', ddim[2], 'samples.')
	} else {
		return(x)
	}

	## subsetting the QC information
	if (!missing(j)) {
		if (!is.null(x@QC)) {
			QC <- x@QC
			if (!is.null(QC$sampleSummary))
				if (ncol(QC$sampleSummary) == ddim[2])
					QC$sampleSummary <- QC$sampleSummary[,j,drop=FALSE]
			if (!is.null(QC$BeadStudioSummary))
				if (nrow(QC$BeadStudioSummary) == ddim[2])
					QC$BeadStudioSummary <- QC$BeadStudioSummary[j,,drop=FALSE]
			x@QC <- QC
		}
		if (!is.null(attr(x, 'vstParameter'))) {
			vstParameter <- attr(x, 'vstParameter')
			if (!is.null(nrow(vstParameter))) {
				if (nrow(vstParameter) == ddim[2]) {
					attr(x, 'vstParameter') <- attr(x, 'vstParameter')[j,,drop=FALSE]
					attr(x, 'transformFun') <- attr(x, 'transformFun')[j]
				}
			}
		}

		## controlData information
		if (nrow(x@controlData) > 0) {
			if (is.numeric(j))  j <- sampleName[j]
			if (all(j %in% colnames(x@controlData)))
				x@controlData <- x@controlData[,j, drop=FALSE]
		}
	}

    # history tracking
    history.finished <- as.character(Sys.time())
	if (is.null(x@history$lumiVersion)) x@history$lumiVersion <- rep(NA, nrow(x@history))
	lumiVersion <- packageDescription('lumi')$Version
	x@history<- rbind(x@history, data.frame(submitted=history.submitted, finished=history.finished, 
			command=history.command, lumiVersion=lumiVersion))

	return(x)
})


setMethod("combine", signature=c(x="LumiBatch", y="LumiBatch"), function(x, y, ...) 
{
	if (missing(y)) return(x)
	if (length(list(...)) > 0) 
	        return(combine(x, combine(y, ...)))
	if (class(x) != class(y))
		stop(paste("objects must be the same class, but are ",
                 class(x), ", ", class(y), sep=""))
	
	if (any(sort(featureNames(x)) != sort(featureNames(y)))) stop('Two data sets have different row names!')
	## determine whether there are duplicated sample names
	sampleName.x <- sampleNames(x)
	sampleName.y <- sampleNames(y)
	if (any(sampleName.x %in% sampleName.y)) {
		warning('Two data sets have some duplicated sample names!\n "_2" were attached to the duplicated sample names!')
		sampleName.x <- sampleNames(x)   # paste(sampleNames(x), '_1', sep='')
		sampleName.y <- paste(sampleNames(y), '_2', sep='')
		sampleNames(x) <- sampleName.x
		sampleNames(y) <- sampleName.y
	}
	featureName.x <- featureNames(x)
	featureName.y <- featureNames(y)
	featureName.com <- featureName.x[featureName.x %in% featureName.y]
	if (length(featureName.com) < length(featureName.x) || length(featureName.com) != length(featureName.y)) {
		if (length(featureName.com) > 0) {
			warning('Two data sets have different featureNames, only the common ones were used!')			
		} else {
			stop('Two data sets have totally different featureNames!')
		}
	}
	## make sure two data sets have the sample order of features
	x <- x[featureName.com,]
	y <- y[featureName.com,]
	history.submitted <- as.character(Sys.time())
	dimm.x <- dim(x); dimm.y <- dim(y)
	assayData(x) <- combine(assayData(x), assayData(y))
	experimentData(x) <- combine(experimentData(x),experimentData(y))

	## combine pheno data
	phenoData(x) <- combine(phenoData(x),phenoData(y))
	
	## featureData(x) <- combine(featureData(x),featureData(y)) # very slow
	## For the feature data, we assume all the data have the same information,
	##    so only the first feature data will be used.
	# if (!is.null(featureData(x)) || !is.null(featureData(y))) {
	# 	feature.x <- featureData(x)
	# 	feature.y <- featureData(y)
    # 
	# 	pData.x <- pData(feature.x)
	# 	pData.y <- pData(feature.y)
	# 	if (names(pData.x)[1] != names(pData.y)[1])	stop('The featureData of two objects are incompatible!')
	# 	#repInfo <- merge(pData.x, pData.y, by=names(pData.x)[1], all=TRUE, suffixes = c(".x",".y"), sort=FALSE)
	# 	#pData(feature.x) <- repInfo
	# 	#
	# 	#metaInfo <- rbind(varMetadata(feature.x), varMetadata(feature.y)[-1,])
	# 	#rownames(metaInfo) <- names(repInfo)
	# 	#varMetadata(feature.x) <- metaInfo
	# 	#featureData(x) <- feature.x
	# }

	## combining the QC information
	if (length(x@QC) > 0 && length(y@QC) > 0) {
		if (!is.null(x@QC$BeadStudioSummary) && !is.null(y@QC$BeadStudioSummary)) {
			if (ncol(x@QC$BeadStudioSummary) == ncol(y@QC$BeadStudioSummary)) {
				BeadStudioSummary <- rbind(x@QC$BeadStudioSummary, y@QC$BeadStudioSummary)
				x@QC$BeadStudioSummary <- BeadStudioSummary
			} else {
				x@QC <- NULL
			}
		} else {
			x@QC <- NULL
		}
		if (!is.null(x@QC$sampleSummary) && !is.null(y@QC$sampleSummary)) {
			if (nrow(x@QC$sampleSummary) == nrow(y@QC$sampleSummary)) {
				sampleSummary <- cbind(x@QC$sampleSummary, y@QC$sampleSummary)
				x@QC$sampleSummary <- sampleSummary
			} else {
				x@QC <- NULL
			}
		} else {
			x@QC <- NULL
		}
		if (!is.null(x@QC)) {
			history.x <- x@QC$history
			if (is.null(history.x)) history.x <- data.frame(submitted=NA, finished=NA, command=NA, lumiVersion=NA)
			if (is.null(history.x$lumiVersion)) history.x$lumiVersion <- rep(NA, nrow(history.x))
			history.y <- y@QC$history
			if (is.null(history.y)) history.y <- data.frame(submitted=NA, finished=NA, command=NA, lumiVersion=NA)
			if (is.null(history.y$lumiVersion)) history.y$lumiVersion <- rep(NA, nrow(history.y))
			x@QC$history <- rbind(history.x, history.y)
		} 
	} else {
		x@QC <- NULL
	}
	
	## VST transformation parameters
	if (!is.null(attr(x, 'vstParameter')) && !is.null(attr(y, 'vstParameter'))) {
		vstParameter.x <- attr(x, 'vstParameter')
		vstParameter.y <- attr(y, 'vstParameter')
		if (is.null(nrow(vstParameter.x))) {
			vstParameter.x <- matrix(vstParameter.x, nrow=1)
		}
		if (is.null(nrow(vstParameter.y))) {
			vstParameter.y <- matrix(vstParameter.y, nrow=1)
		}
		if (nrow(vstParameter.x) != dimm.x[2] || nrow(vstParameter.y) != dimm.y[2]) {
			attr(x, 'vstParameter') <- attr(x, 'transformFun') <- NULL
		} else {
			attr(x, 'vstParameter') <- rbind(attr(x, 'vstParameter'), attr(y, 'vstParameter'))
			attr(x, 'transformFun') <- c(attr(x, 'transformFun'), attr(y, 'transformFun'))
		}
	}
	
	## controlData information
	if (nrow(x@controlData) > 0) {
		if (nrow(x@controlData) == nrow(y@controlData)) {
			controlData <- cbind(x@controlData, y@controlData)
			x@controlData <- as.data.frame(controlData)
		}
	}

	# history tracking
	history.finished <- as.character(Sys.time())
	history.command <- capture.output(print(match.call(combine)))  
	x@history<- rbind(x@history, y@history)
	if (is.null(x@history$lumiVersion)) x@history$lumiVersion <- rep(NA, nrow(x@history))
	lumiVersion <- packageDescription('lumi')$Version
    x@history<- rbind(x@history, 
	       data.frame(submitted=history.submitted, finished=history.finished, command=history.command, lumiVersion=lumiVersion))
	return(x)
})


##some special handling of main is needed
setMethod("boxplot",signature(x="ExpressionSet"),
	function(x, range=0, main, logMode=TRUE, subset=5000, ...) 
{
  	tmp <- description(x)
  	if (missing(main) && (is(tmp, "MIAME")))
     	main <- tmp@title

	expr <- exprs(x)
  	if (logMode && max(exprs(x), na.rm=TRUE) > 50) {
		## force the expression value as positive in the logMode
		# if (min(expr, na.rm=TRUE) < 0) expr <- expr - min(expr, na.rm=TRUE) + 1
		# remove the negative values
		if (min(expr) < 0) {
			rMin <- rowMin(expr)
			expr <- expr[rMin > 0, , drop=FALSE]
		}
		expr <- log2(expr)
	} 
	if (!is.null(subset)) {
		if (!is.numeric(subset)) stop('subset should be numeric!')
		if (length(subset) == 1) {
			index <- sample(1:nrow(expr), min(subset, nrow(expr)))
		} else {
			index <- subset
			index <- index[index > 0 & index <= nrow(expr)]
		}
	} else {
		index <- 1:nrow(expr)
	}

	dataMatrix <- expr[index,]
	labels <- colnames(dataMatrix)
	if (is.null(labels)) labels <- as.character(1:ncol(dataMatrix))
	## set the margin of the plot
	mar <- c(max(nchar(labels))/2 + 4.5, 5, 5, 3)
	old.mar <- par('mar')
	old.xaxt <- par('xaxt')
	par(xaxt='n')
	par(mar=mar)
	boxplot(dataMatrix ~ col(dataMatrix), main=main, range=range, xlab='', ylab='amplitude', ...)
	par(xaxt='s')
	axis(1, at=1:ncol(dataMatrix), labels=labels, tick=TRUE, las=2)
	par(mar=old.mar)
	par(xaxt=old.xaxt)
})


setMethod('hist', signature(x='ExpressionSet'), 
	function(x, ...) 
{
	density(x, ...)
})


setMethod('density', signature(x='ExpressionSet'), 
	function(x, logMode=TRUE, xlab = NULL, ylab = "density", type = "l", col=1:dim(x)[2], lty=1:dim(x)[2], 
	lwd=1, xlim=NULL, index.highlight=NULL, color.highlight=2, symmetry=NULL, addLegend=TRUE, subset=5000, main='',...) 
{
	if (is(x, 'ExpressionSet')) {
	    expr <- exprs(x)
	} else if (is.numeric(x)) {
		expr <- as.matrix(x)
	} else {
		stop('Un-supported class of x!')
	}
		
    if (logMode && (max(expr, na.rm=TRUE) > 50)) {
		## force the expression value as positive in the logMode
		# if (min(expr, na.rm=TRUE) < 0) expr <- expr - min(expr, na.rm=TRUE) + 1
		# remove the negative values
		if (min(expr) < 0) {
			rMin <- rowMin(expr)
			expr <- expr[rMin > 0, , drop=FALSE]
		}
		expr <- log2(expr)
		if (is.null(xlab)) 
			xlab <- "log2 intensity"
    } else if (is.null(xlab)) 
        xlab <- "intensity"

	if (!is.null(subset)) {
		if (!is.numeric(subset)) stop('subset should be numeric!')
		if (length(subset) == 1) {
			index <- sample(1:nrow(expr), min(subset, nrow(expr)))
		} else {
			index <- subset
			index <- index[index > 0 & index <= nrow(expr)]
		}
	} else {
		index <- 1:nrow(expr)
	}
	expr <- expr[index,,drop=FALSE]

	if (!is.null(symmetry)) {
		x.range <- range(expr)
		if (symmetry > x.range[1] && symmetry < x.range[2]) {
			warning('symmetry point should not be within the range of x!')
			symmetry <- NULL
		} else {
			expr <- rbind(expr, 2*symmetry - expr)
		}
	}
	x.density <- apply(expr, 2, density, ...)
    all.x <- do.call("cbind", lapply(x.density, function(x) x$x))
    all.y <- do.call("cbind", lapply(x.density, function(x) x$y))

	if (!is.null(symmetry)) {
		nr <- nrow(all.x)
		if (all.x[1,1] >= x.range[1] && all.x[1,1] <= x.range[2]) {
			all.x <- all.x[1:round(nr/2),]
			all.y <- all.y[1:round(nr/2),]
		} else {
			all.x <- all.x[round(nr/2):nr,]
			all.y <- all.y[round(nr/2):nr,]
		}
	}
	# matplot(all.x, all.y, ylab=ylab, xlab=xlab, type=type, col=col, lty=lty, lwd=lwd, ...)
	matplot(all.x, all.y, ylab=ylab, xlab=xlab, type=type, col=col, lty=lty, lwd=lwd, xlim=xlim, main=main)
	if (!is.null(index.highlight)) {
		if (index.highlight > ncol(all.x) || index.highlight < 1) {
			warning('Highlight index out of range!')
			index.highlight <- 1
		}
		lines(all.x[,index.highlight], all.y[,index.highlight], col=color.highlight, lwd=2, lty=1)
	}
	## add legend
	if (addLegend) {
		labels <- colnames(expr)
		if (is.null(labels)) labels <- as.character(1:ncol(expr))

		col <- col
		lwd <- lwd
		lty <- lty
		if (!is.null(index.highlight)) {
			col[index.highlight] <- color.highlight
			lwd[index.highlight] <- 2
			lty[index.highlight] <- 1
		}
		x.pos <- (max(all.x) - min(all.x)) * 2/3 + min(all.x)
		y.pos <- max(all.y)
		legend(x.pos, y.pos, legend=labels, col=col, lwd=lwd, lty=lty)
	}
	
    invisible(list(all.x = all.x, all.y = all.y))
})


setMethod("pairs", signature(x="ExpressionSet"), 
	function(x,..., logMode=TRUE, subset=5000) 
{
	upperPanel <- function(x, y, fold=2) {
		points(x[subset], y[subset])
		abline(0, 1, col="red", lty=1)
		if (logMode) {
			abline(log2(fold), 1, col="green", lty=2)
			abline(log2(1/fold), 1, col="green", lty=2)
		} else {
			abline(fold, 1, col="green", lty=2)
			abline(-fold, 1, col="green", lty=2)
		}
	}

	lowerPanel <- function(x, y, cex=1.44, fold=2) {
		if (logMode) {
			up <- length(which((x-y) > log2(fold)))
			down <- length(which((y-x) > log2(fold)))
		} else {
			up <- length(which((x/y) > fold))
			down <- length(which((y/x) > fold))
		}
		ex <- par("fin")[1]*0.9
		txt <- paste("Cor =", as.character(round(cor(x,y),2)),"\n")
		txt <- paste(txt, up, " (> ", fold, ", up)\n", sep="")
		txt <- paste(txt, down, " (> ", fold, ", down)\n", sep="")
		text(mean(range(x)), mean(range(x)), labels=txt, cex=ex)
	}

	## put histograms on the diagonal
	diagPanel <- function(x, ...) {
	    usr <- par("usr"); on.exit(par(usr))
	    par(usr = c(usr[1:2], 0, 1.5) )
	    h <- hist(x, plot = FALSE)
	    breaks <- h$breaks; nB <- length(breaks)
	    y <- h$counts; y <- y/max(y)
	    rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
	}

	expr <- exprs(x)
	if(logMode) {
		if (max(expr, na.rm=TRUE) > 50) {
			## force the expression value as positive in the logMode
			# if (min(expr, na.rm=TRUE) < 0) expr <- expr - min(expr, na.rm=TRUE) + 1
			# remove the negative values
			if (min(expr) < 0) {
				rMin <- rowMin(expr)
				expr <- expr[rMin > 0, , drop=FALSE]
			}
			expr <- log2(expr)
		}
	} else {
		if (max(expr, na.rm=TRUE) < 50) {
			expr <- 2^expr
		}
	}

	if (!is.null(subset)) {
		if (!is.numeric(subset)) stop('subset should be numeric!')
		if (length(subset) == 1) {
			subset <- sample(1:nrow(expr), min(subset, nrow(expr)))
		} else {
			subset <- subset[subset > 0 & subset <= nrow(expr)]
		}
	} else {
		subset <- 1:nrow(expr)
	}
	
	pairs(expr,upper.panel=upperPanel, diag.panel=diagPanel, 
			lower.panel=lowerPanel, ...)
})


if(is.null(getGeneric("MAplot")))
	setGeneric("MAplot", function(object, ...)
		standardGeneric("MAplot"))
  	

setMethod("MAplot", signature(object="ExpressionSet"), 
	function(object, ..., logMode=TRUE, subset=5000) 
{
	expr <- exprs(object)
	if(logMode) {
		if (max(expr, na.rm=TRUE) > 50) {
			## force the expression value as positive in the logMode
			# if (min(expr, na.rm=TRUE) < 0) expr <- expr - min(expr, na.rm=TRUE) + 1
			# remove the negative values
			if (min(expr) < 0) {
				rMin <- rowMin(expr)
				expr <- expr[rMin > 0, ,drop=FALSE]
			}
			expr <- log2(expr)
		} 
	} else {
		if (max(expr, na.rm=TRUE) < 50) {
			expr <- 2^expr
		}
	}
	if (!is.null(subset)) {
		if (!is.numeric(subset)) stop('subset should be numeric!')
		if (length(subset) == 1) {
			index <- sample(1:nrow(expr), min(subset, nrow(expr)))
		} else {
			index <- subset
			index <- index[index > 0 & index <= nrow(expr)]
		}
	} else {
		index <- 1:nrow(expr)
	}
	mva.pairs(expr[index, ], log.it=FALSE, ...)
})


setMethod('plot',
	signature('LumiBatch', 'missing'),
	function(x, what=c('density', 'boxplot', 'pair', 'MAplot', 'sampleRelation', 'outlier', 'cv'), main, ...)
{

	object <- x
	if (!is(object, 'LumiBatch')) stop('The object should be class "LumiBatch"!')
	what <- match.arg(what)
	
	if (what == 'density') {
		if (missing(main)) main <- 'Density plot of intensity'
		density(object, xlab="intensity", ylab="density", main=main, ...)
	} else if (what == 'boxplot') {
		if (missing(main)) main <- 'Boxplot of microarray intensity'
		boxplot(object, xlab='microarrays', ylab='intensity', main=main, ...)
	} else if (what == 'cv') {
		if (missing(main)) main <- 'Density plot of coefficient of variance'
		estimateLumiCV(object, ifPlot=TRUE, main=main, ...)
	} else if (what == 'sampleRelation') {
		plotSampleRelation(object, ...)
	} else if (what == 'pair') {
		if (missing(main)) main <- 'Pairwise plot with sample correlation'
		pairs(object, main=main, ...)
	} else if (what == 'MAplot') {
		if (missing(main)) main <- 'Pairwise MA plots between samples'
		MAplot(object, main=main, ...)
	} else if (what == 'outlier') {
		detectOutlier(object, ifPlot=TRUE, ...)
	} else {
		print('Unsupported !')
	}
	return(invisible(TRUE))	
})
