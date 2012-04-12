lumiR.batch <- function(fileList, convertNuID = TRUE, lib = NULL, transform = c('none', 'vst', 'log2', 'cubicRoot'), sampleInfoFile = NULL, ...) {

	oldDir <- getwd()
	dirMode <- FALSE
	transform <- match.arg(transform)
	if (file.exists(fileList[1])) {
		if (file.info(fileList[1])[1,'isdir']) {
			dirMode <- TRUE
			setwd(fileList)
		}
	}
	if (dirMode && length(fileList) == 1) {
		fileList <- dir(fileList, pattern='.csv')
		if (length(fileList) == 0) stop('No data files were found!')
	}

	print('Inputting the data ...')
	for (i in 1:length(fileList)) {
		file.i <- fileList[i]
		x.lumi.i <- lumiR(file.i, parseColumnName=FALSE, convertNuID = FALSE, ...)
		# x.lumi.i <- lumiR(file.i, parseColumnName=FALSE)
		if (transform != 'none') {
			x.lumi.i <- lumiT(x.lumi.i, method=transform, simpleOutput=TRUE)
		}
		if (i == 1) {
			x.lumi <- x.lumi.i
		} else {
			x.lumi <- combine(x.lumi, x.lumi.i)
		}
	}
	if (!convertNuID) lib <- NULL	
	if (!is.null(lib) || convertNuID) {
		print('Adding nuID to the data ...')
		x.lumi <- addNuId2lumi(x.lumi, lib=lib)
	}

	if (!is.null(sampleInfoFile)) {

		if (is.character(sampleInfoFile) || class(sampleInfoFile)[1] == 'file') {
			if (file.exists(sampleInfoFile)) {
				sampleInfo <- read.table(sampleInfoFile, head=TRUE, sep='\t', colClasses='character', comment='')
			} else {
				warning('The provided sampleInfoFile does not exist!')
				setwd(oldDir)
				return(x.lumi)
			}
		} else if (is.data.frame(sampleInfoFile)) {
			sampleInfo <- sampleInfoFile
		}
		## force the names to be capitalized
		names(sampleInfo) <- toupper(names(sampleInfo))
		ID <- sampleInfo$ID
		if (is.null(ID)) {
			ID <- sampleInfo[,1]
			if (any(duplicated(ID))) {
				warning('In sampleInfoFile, the ID column is required or the first column should be unique!')
				setwd(oldDir)
				return(x.lumi)
			}
		} 
		rownames(sampleInfo) <- ID

		sampleName <- sampleNames(x.lumi)
		ID <- ID[ID %in% sampleName]
		if (nrow(sampleInfo) != length(ID)) {
			warning('Some IDs provided in the sampleInfoFile do not exist the data file!')
			if (length(ID) == 0) {
				stop('The IDs provided in the sampleInfoFile do not match the data file!')
			} 
			
		} 
		x.lumi <- x.lumi[, ID]

		pData <- sampleInfo[ID,]
		label <- sampleInfo[ID, 'LABEL']
		if (!is.null(label)) {
			label <- ID
			rownames(pData) <- label
			sampleNames(x.lumi) <- label
		}

		pdata <- new("AnnotatedDataFrame", data=pData)
		phenoData(x.lumi) <- pdata
	} 
	setwd(oldDir)
	return(x.lumi)
}