`lumiR` <-
function(fileName, sep = NULL, detectionTh = 0.01, na.rm = TRUE, convertNuID = TRUE, lib = NULL, dec='.', parseColumnName=FALSE, checkDupId=TRUE, 
	columnNameGrepPattern=list(exprs='AVG_SIGNAL', se.exprs='BEAD_STD', detection='Detection', beadNum='Avg_NBEADS'),
	inputAnnotation=TRUE, annotationColumn=c('ACCESSION', 'SYMBOL', 'PROBE_SEQUENCE', 'PROBE_START', 'CHROMOSOME', 'PROBE_CHR_ORIENTATION', 'PROBE_COORDINATES', 'DEFINITION'), ...) 
{
	## the patterns used to grep columns in the BeadStudio output text file 
	## 'exprs' and 'se.exprs' related columns are required
	if (is.null(columnNameGrepPattern$exprs)) columnNameGrepPattern$exprs <- 'AVG_SIGNAL'
	if (is.null(columnNameGrepPattern$se.exprs)) columnNameGrepPattern$se.exprs <- 'BEAD_STD'
	if (is.null(columnNameGrepPattern$detection)) columnNameGrepPattern$detection <- 'Detection'
	if (is.null(columnNameGrepPattern$beadNum)) columnNameGrepPattern$beadNum <- 'Avg_NBEADS'

	if (is.na(columnNameGrepPattern$exprs)) {
		columnNameGrepPattern$exprs <- 'AVG_SIGNAL'
		warning('exprs slot is required and default pattern will be used!')
	}
	if (is.na(columnNameGrepPattern$se.exprs)) {
		columnNameGrepPattern$beadNum <- columnNameGrepPattern$detection <- NA
		warning('se.exprs slot is required for the VST transformation!\n We strongly suggest to include BEAD_STD columns!')
		# columnNameGrepPattern$se.exprs <- 'BEAD_STD'
		# warning('se.exprs slot is required and default pattern will be used!')
	}
	if (is.na(columnNameGrepPattern$beadNum)) {
		# columnNameGrepPattern$beadNum <- 'Avg_NBEADS'
	}

	# columnNameGrepPattern <- c(exprs='AVG_SIGNAL', se.exprs='BEAD_STD', detection='Detection', beadNum='Avg_NBEADS')
	history.submitted <- as.character(Sys.time())
	## set "stringsAsFactors" as FALSE
	oldSetting <- options("stringsAsFactors")[[1]]
	options(stringsAsFactors = FALSE)
	version <- 2

	if (!file.exists(fileName)) stop('The file is not exist! Please check your file path!')
	## ---------------------------------------
	## identify the Metadata lines 
	info <- readLines(fileName, n=20)    # take the first 20 lines to have a taste

	## Use "AVG_SIGNAL" as an indicator of Where the metaData stops
	##   intelligently find nMetaDataLines  
	nMetaDataLines <- grep(columnNameGrepPattern$exprs, info, ignore.case=TRUE) - 1
    
	if (is.null(sep)) {
	    ## Find out the separator (sep) by taking the first two line of data, and comparing them.
	    ##  we assume it is either "," or "\t".
    	titleLine <- info[nMetaDataLines + 1]
	    dataLine1 <- info[nMetaDataLines + 2]
		dataLine2 <- info[nMetaDataLines + 3]
		sepNum1 <- gregexpr('\t', dataLine1)[[1]]
		sepNum2 <- gregexpr('\t', dataLine2)[[1]]
		if (sepNum1[1] > 0 && length(sepNum1) == length(sepNum2)) {
			sep <- '\t'
		} else if (dec != ',') {
			sepNum1 <- gregexpr(',', dataLine1)[[1]]
			sepNum2 <- gregexpr(',', dataLine2)[[1]]
			if (sepNum1[1] > 0 && length(sepNum1) == length(sepNum2)) {
				sep <- ','
			} else {
				stop('The seperator is not Tab or comma!\n Please sepecify the seperator used in the file!')
			}
		} else {
			stop('Please sepecify the seperator used in the file!')
		}
	}
	## determine whether the quote is used or not
	dataLine1 <- strsplit(info[nMetaDataLines + 2], sep)[[1]]
	quoteCount1 <- gregexpr('"', dataLine1[1])[[1]]
	quoteCount2 <- gregexpr('\'', dataLine1[1])[[1]]
	if (length(quoteCount1) == 2) {
		quote <- '"'
	} else if (length(quoteCount2) == 2) {
		quote <- '\''
	} else {
		quote <- ''
	}

	header <- strsplit(info[nMetaDataLines + 1], sep)[[1]]
	probeId.pos <- grep('Probe.?ID', header, ignore.case=TRUE)
	if (length(probeId.pos) > 0) {
		colClasses <- rep(NA, length(header))
		colClasses[probeId.pos] <- 'character'
	} else {
		colClasses <- NULL
	}

	## ---------------------------------------
	# get meta data info
	if (nMetaDataLines > 0) {
		info <- readLines(fileName, n=nMetaDataLines)
		## check the version of the beadStudio output
		markerInd <- grep('^\\[.*\\]', info, ignore.case=TRUE)
		if (length(markerInd) > 0) {
			if (length(grep('^\\[Header\\]', info[markerInd[1]], ignore.case=TRUE)) == 0) 
				warning('The data file may not be in the Illumia BeadStudio output format!')
			if (length(markerInd) > 1) {
				if (length(grep('^\\[.*\\Profile]', info[markerInd[2]], ignore.case=TRUE)) == 0) 
					warning('The data file may not be in the Illumia BeadStudio output format!')
			}
			version <- 3
			info <- info[-markerInd]
		}
		if (length(info) > 0) {
			## remove the blanks
			info <- sub("[[:blank:]]+$", "", info)
			info <- sub(paste(sep,"+$", sep=''), "", info)

			## check the meta info of the file
			if (version == 2) {
				ind <- grep("BeadStudio version", info, ignore.case=TRUE)
			} else {
				ind <- grep("BSGX Version", info, ignore.case=TRUE)
			}
			if (length(ind) == 0) 
			    warning("The data file may not be in the Illumia BeadStudio output format.")

			## should not be normalized in BeadStudio
			ind <- grep("Normalization", info, ignore.case=TRUE)  # find where is the row index
			if (length(ind) > 0) {
				if (version == 2) {
					normalization <- strsplit(info, split='=')[[ind]][2]
					normalization <- gsub(pattern=" |,", replace="", normalization) # remove space or ","
				} else {
					normalization <- strsplit(info, split=sep)[[ind]][2]
				}
				if (length(grep("none", normalization, ignore.case=TRUE)) == 0) {
				    warning("The raw data should not be normalized in BeadStudio.")
				}
			}
		} else {
			info <- NULL
		}
	} else {
		info <- NULL
	}
    
	allData <- read.table(file=fileName, header=TRUE, sep=sep, dec = dec, skip=nMetaDataLines, row.names=NULL, colClasses=colClasses,
		quote=quote, as.is=TRUE, check.names=FALSE, strip.white=TRUE, comment.char="", fill=TRUE, ...)

	## retrieve the possible section line index
	sectionInd <- grep('^\\[.*\\]', allData[,1], ignore.case=TRUE)
    
	if (length(sectionInd) > 0) {
		if (is.na(version)) verion <- 3
		otherData <- allData[sectionInd[1]:nrow(allData), ]
		## we assume the first section is the expression data section
		allData <- allData[1:(sectionInd[1]-1),, drop=FALSE]
		## remove the all NA columns, which can be produced when saved in Excel
		naCol <- apply(allData, 2, function(x) all(is.na(x) | x == ''))
		allData <- allData[,!naCol]
    
		## process otherData
		sectionInd <- sectionInd - sectionInd[1] + 1
		sectionName <- otherData[sectionInd, 1]
    
		## retrieve the control data
		controlInd <- grep('^\\[Control.*\\]', sectionName, ignore.case=TRUE)
		if (length(controlInd) > 1) {
			ind <- grep('^\\[Control probe.*\\]', sectionName[controlInd], ignore.case=TRUE)
			if (length(ind) > 0) {
				controlInd <- controlInd[ind[1]]
			} else {
				controlInd <- controlInd[1]
			}
		}
		if (length(controlInd) > 0) {
			startRow <- sectionInd[controlInd] + 1
			if (length(sectionInd) > controlInd)
				endRow <- sectionInd[controlInd + 1] - 1
			else 
				endRow <- nrow(otherData)
			controlData <- otherData[startRow:endRow,]
			## remove the all NA columns, which can be produced when save in Excel
			naCol <- apply(controlData, 2, function(x) all(is.na(x) | x == ''))
			controlData <- controlData[,!naCol, drop=FALSE]
			colnames(controlData) <- controlData[1,]
			controlData <- controlData[-1,]
		} else {
			controlData <- data.frame()
		}
    
		## retrieve the Sample Table
		summaryInd <- grep('^\\[Sample.*Table\\]', sectionName, ignore.case=TRUE)
		if (length(summaryInd) > 0) {
			startRow <- sectionInd[summaryInd] + 1
			if (length(sectionInd) > summaryInd)
				endRow <- sectionInd[summaryInd + 1] - 1
			else 
				endRow <- nrow(otherData)
			sampleSummary <- otherData[startRow:endRow,]
			## remove the all NA columns, which can be produced when save in Excel
			naCol <- apply(sampleSummary, 2, function(x) all(is.na(x) | x == ''))
			sampleSummary <- sampleSummary[,!naCol]
			colnames(sampleSummary) <- sampleSummary[1,]
			sampleSummary <- sampleSummary[-1,]
		} else {
			sampleSummary <- data.frame()
		}
	} else {
		controlData <- sampleSummary <- data.frame()
	}
	header <- names(allData)
    
	## Get Id. The ProbeID (by default it is the second column) is preferred if provided, 
	# otherwise the TargetID (by default it is the first column) is used.
	targetID <- as.character(as.vector(allData[,1]))
	id <- targetID
	idName <- header[1]
	if (length(grep('Probe.?ID', header[2], ignore.case=TRUE)) > 0) {
		id <- as.character(as.vector(allData[,2]))
		idName <- header[2]
	} else if (!is.null(lib)) {
		probeId.pos <- grep('Probe.?ID', header, ignore.case=TRUE)
		if (length(probeId.pos) > 0) {
			id <- as.character(as.vector(allData[,probeId.pos[1]]))
			idName <- header[probeId.pos[1]]
		}
	}
    
	## identify where the signal column exists
	ind <- grep(columnNameGrepPattern$exprs, header, ignore.case=TRUE)
	if (length(ind) == 0) stop('Input data format unrecognizable!\nThere is no column name contains "AVG_SIGNAL"!')
	exprs <- as.matrix(allData[,ind])
	if (!is.numeric(exprs[1])) {
		exprs <- matrix(as.numeric(exprs), nrow=nrow(allData))
	} 
	colnames(exprs) <- header[ind]
	
	## identify where the signal standard deviation column exists 
	if (is.na(columnNameGrepPattern$'se.exprs')) {
		ind <- NULL
	} else {
		ind <- grep(columnNameGrepPattern$'se.exprs', header, ignore.case=TRUE)
	}
	if (length(ind) == 0) {
		se.exprs <- NULL
		 # stop('Input data format unrecognizable!\nThere is no column name contains "BEAD_STDEV"!')
	} else {
		se.exprs <- as.matrix(allData[,ind])
		if (!is.numeric(se.exprs[1])) {
			se.exprs <- matrix(as.numeric(se.exprs), nrow=nrow(allData))
		}
		colnames(se.exprs) <- header[ind]
	}
	## identify the detection columns
	if (is.na(columnNameGrepPattern$detection)) {
		ind <- NULL
	} else {
		ind <- grep(columnNameGrepPattern$detection, header, ignore.case=TRUE)
	}
	if (length(ind) == 0) {
		detection <- NULL
	} else {
		detection <- as.matrix(allData[,ind])
		if (!is.numeric(detection[1])) {
			detection <- matrix(as.numeric(detection), nrow=nrow(allData))
		}
		colnames(detection) <- header[ind]
    
		if (length(grep("Detection Pval", header, ignore.case=TRUE)) == 0) {
			detection <- 1 - detection
		}
	}
	## identify the bead number columns
	if (is.na(columnNameGrepPattern$beadNum)) {
		ind <- NULL
	} else {
		ind <- grep(columnNameGrepPattern$beadNum, header, ignore.case=TRUE)
	}
	if (length(ind) == 0) {
		beadNum <- NULL
	} else {
	    beadNum <- as.matrix(allData[,ind])
		if (!is.numeric(beadNum[1])) {
			beadNum <- matrix(as.numeric(beadNum), nrow=nrow(allData))
		} 
		colnames(beadNum) <- header[ind]
	}
    
	## identify the annotation columns
	annotationInfo <- NULL
	if (inputAnnotation) {
		# It is based on annotationColumn
		annotationColumn <- header[toupper(header) %in% toupper(annotationColumn)]
		if (length(annotationColumn) == 0 && "SPECIES" %in% toupper(header)) {
			print('Some annotation columns not available in the data.')
		} else {
			annotationInfo <- allData[,annotationColumn]
		}
	}

	## check for possible duplicated ids
	if (checkDupId) {
		dupId <- unique(id[duplicated(id)])
		if (length(dupId) > 0) {
			warning('Duplicated IDs found and were merged!')
			rmInd <- NULL
			for (dupId.i in dupId) {
				selInd.i <- which(id == dupId.i)
				exprs[selInd.i[1],] <- colMeans(exprs[selInd.i,])
				if (is.null(beadNum)) {
					se.exprs[selInd.i[1],] <- colMeans(se.exprs[selInd.i,])
				} else {
					totalBead.i <- colSums(beadNum[selInd.i,])
					beadNum[selInd.i[1],] <- totalBead.i				
					temp <- colSums(se.exprs[selInd.i,]^2 * (beadNum[selInd.i,] - 1))
					temp <- temp / (totalBead.i - length(selInd.i))
					se.exprs[selInd.i[1],] <- sqrt(temp * (colSums(1/beadNum[selInd.i,])))
				}
				if (!is.null(detection)) {
					detection[selInd.i[1],] <- apply(detection[selInd.i,], 2, max)
				}
				rmInd <- c(rmInd, selInd.i[-1])
			}
			## remove duplicated
			exprs <- exprs[-rmInd,,drop=FALSE]
			if (!is.null(se.exprs)) se.exprs <- se.exprs[-rmInd,,drop=FALSE]
			id <- id[-rmInd]; targetID <- targetID[-rmInd]
			if (!is.null(detection)) detection <- detection[-rmInd,,drop=FALSE]
			if (!is.null(beadNum)) beadNum <- beadNum[-rmInd,,drop=FALSE]
			if (!is.null(annotationInfo)) annotationInfo <- annotationInfo[-rmInd,,drop=FALSE]
		}
	}
    
	if (na.rm) {
		## remove the probes with all of them as NA
		keepInd <- apply(is.na(exprs), 1, sum) == 0
		exprs <- exprs[keepInd,,drop=FALSE]
		se.exprs <- se.exprs[keepInd,,drop=FALSE]
		if (!is.null(beadNum))  beadNum <- beadNum[keepInd,,drop=FALSE]
		if (!is.null(detection))  detection <- detection[keepInd,,drop=FALSE]
		if (!is.null(annotationInfo)) annotationInfo <- annotationInfo[keepInd,,drop=FALSE]
		id <- id[keepInd]
		targetID <- targetID[keepInd]
	}
	if (checkDupId) {
		rownames(exprs) <- id
		if (!is.null(se.exprs)) rownames(se.exprs) <- id
		if (!is.null(beadNum)) rownames(beadNum) <- id
		if (!is.null(detection)) rownames(detection) <- id
	}
    
	# get sample information
	pattern <- paste('[^[:alnum:]]*', columnNameGrepPattern$exprs, '[^[:alnum:]]*', sep='')
	sampleID <-  sub(pattern, '', colnames(exprs), ignore.case=TRUE) 
	if (any(duplicated(sampleID))) {
		warning('Duplicated column names found in the raw data! \n A suffix number is added to the duplicated column names.')
		dupId <- which(duplicated(sampleID))
		dupName <- unique(sampleID[dupId])
		for (dupName.i in dupName) {
			dupInd.i <- which(sampleID == dupName.i)
			sampleID[dupInd.i] <- paste(sampleID[dupInd.i], 1:length(dupInd.i), sep='.')
		}
	}
	if (parseColumnName) {
		sampleIDInfo <- strsplit(sampleID, split="_")
		label <- NULL
		newID <- NULL
		temp <- lapply(sampleIDInfo, function(x) {
			label <<- c(label, x[length(x)])
			newID <<- c(newID, paste(x[1:2], collapse="_"))
			})
		if (!any(duplicated(newID))) sampleID <- newID
		if (length(unique(label)) != length(label) || length(label) == 0 || any(is.na(label)))
			label <- sampleID
	} else {
		label <- sampleID
	}
    
	## reportInfo save the id information
	reporterInfo <- data.frame(id)
	names(reporterInfo) <- idName
	varMetadata <- data.frame(labelDescription='The Illumina microarray identifier')
	varName <- idName

	## if ProbeID is used as id, then also keep the TargetID information in the featureData
	if (idName != header[1]) {
		reporterInfo <- data.frame(reporterInfo, TargetID=targetID)
		varMetadata <- data.frame(labelDescription=c(varMetadata$labelDescription, 'The Illumina TargetID'))
		varName <- c(varName, 'TargetID')
	}
	## add annotationInfo to the featureData
	if (!is.null(annotationInfo)) {
		reporterInfo <- data.frame(reporterInfo, annotationInfo)
		varMetadata <- data.frame(labelDescription=c(varMetadata$labelDescription, names(annotationInfo)))
		varName <- c(varName, names(annotationInfo))
	}
	rownames(varMetadata) <- varName
	if (checkDupId)	rownames(reporterInfo) <- id
	featureData <- new("AnnotatedDataFrame", data=reporterInfo, varMetadata=varMetadata)
    
	## check the dimensions of the input data
	colnames(exprs) <- label
	if (!is.null(se.exprs)) {
		if (ncol(exprs) != ncol(se.exprs)) 
			stop('Different column numbers of exprs and se.exprs! Please check the input data format.')
		colnames(se.exprs) <- label
	}
	if (!is.null(beadNum)) {
		if (ncol(beadNum) == length(label)) {
			colnames(beadNum) <- label
		} else {
			warning('The number of beadNum columns does not match! Please check the input data format.')
			if (ncol(beadNum) > length(label)) beadNum <- beadNum[,1:length(label)]
			if (ncol(beadNum) < length(label)) {
				for (i in 1:(length(label) - ncol(beadNum)))
					beadNum <- cbind(beadNum, rep(NA, nrow(beadNum)))
			}
		}
	}
	if (!is.null(detection)) {
		if (ncol(detection) == length(label)) {
			colnames(detection) <- label
		} else {
			warning('The number of detection columns does not match! Please check the input data format.')
			if (ncol(detection) > length(label)) beadNum <- beadNum[,1:length(label)]
			if (ncol(detection) < length(label)) {
				for (i in 1:(length(label) - ncol(detection)))
					detection <- cbind(detection, rep(NA, nrow(detection)))
			}
		}
	}
	## If no se.exprs imported, it will create a ExpressionSet class, instead of LumiBatch class.
	if (is.null(se.exprs)) {
		cmd <- 'x.lumi <- new("ExpressionSet", exprs=exprs'
	} else {
		cmd <- 'x.lumi <- new("LumiBatch", exprs=exprs, se.exprs=se.exprs'
		if (!is.null(detection)) cmd <- paste(cmd, ', detection=detection')
		if (!is.null(beadNum)) cmd <- paste(cmd, ', beadNum=beadNum')
		cmd <- paste(cmd, ', featureData=featureData')		
	}
	
	## produce the phenoData object
	if (!all(sampleID == label)) {
		pData <- data.frame(sampleID=sampleID, label=label)
		rownames(pData) <- label
		varMetadata <- data.frame(labelDescription=c('The unique Illumina microarray Id', 
			'The label of the sample'))
		rownames(varMetadata) <- c('sampleID', 'label')
	}  else {
		pData <- data.frame(sampleID=sampleID)
		rownames(pData) <- sampleID
		varMetadata <- data.frame(labelDescription=c('The unique Illumina microarray Id'))
		rownames(varMetadata) <- c('sampleID')
	}
	pdata <- new("AnnotatedDataFrame", data=pData, varMetadata=varMetadata)
	cmd <- paste(cmd, ', phenoData=pdata')
	cmd <- paste(cmd, ')')
	eval(parse(text=cmd))
	if (is.null(se.exprs)) {
		## resume the old settings
		options(stringsAsFactors = oldSetting)
		return(x.lumi)
	}
	x.lumi@controlData <- controlData
	x.lumi@QC <- list(BeadStudioSummary=sampleSummary)
	sampleNames(x.lumi) <- label
	
	if (!is.null(info)) {
		info <- gsub('\t+', '\t', info)
	}
	experimentData(x.lumi)@other <- list(info)

    # history tracking
    history.finished <- as.character(Sys.time())
	#history.command <- match.call()
    history.command <- capture.output(print(match.call(lumiR)))  
	
	## replace with the real file name
	if (length(grep(',', history.command)) > 0) {
		history.command <- sub('\\(.+,', paste('("', fileName, '",', sep=''), history.command)
	} else {
		history.command <- sub('\\(.+\\)', paste('("', fileName, '")', sep=''), history.command)
	}

	lumiVersion <- packageDescription('lumi')$Version
	x.lumi@history<- data.frame(submitted=history.submitted, 
			finished=history.finished, command=history.command, lumiVersion=lumiVersion)

	## Add the species information if exists
	if (any(toupper(header) == 'SPECIES')) {
		species <- as.character(allData[1,header[toupper(header) == 'SPECIES']])
		annotation(x.lumi) <- species
	}
	
	## initialize the QC slot in the LumiBatch object
	x.lumi <- lumiQ(x.lumi, detectionTh=detectionTh)

	## Add nuID if the annotation library is provided
	if (!convertNuID) lib <- NULL
	if (!is.null(lib) || convertNuID)  x.lumi <- addNuId2lumi(x.lumi, lib=lib)

	## resume the old settings
	options(stringsAsFactors = oldSetting)
    
    return(x.lumi)
}

