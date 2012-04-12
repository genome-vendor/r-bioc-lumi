`addNuId2lumi` <-
function(x.lumi, annotationFile=NULL, sep=NULL, lib=NULL, annotationColName=c(sequence='Probe_Sequence', target='Target', probe='ProbeId')) {

    history.submitted <- as.character(Sys.time())

	## check whether the object is nuID annotated.
	exprs <- exprs(x.lumi)
	id <- rownames(exprs)
	if(all(sapply(id[1:20], is.nuID))) {
		print('The lumiBatch object is already nuID annotated!')
		return(x.lumi)
	}
	if (!is.null(lib)) {
		if (length(grep('\\.db', lib)) > 0) {
			warning(paste(lib, 'does not include nuID conversion information!'))
			lib <- NULL
		}
	}

	newId <- id
	## ---------------------------------------
	## identify the Metadata lines 
	if (!is.null(annotationFile)) {
		info <- readLines(file(annotationFile), n=10)    # take the first 10 lines to have a taste

		## Use annotationColName[1] as an indicator of Where the metaData stops
		##   intelligently find nMetaDataLines  
		nMetaDataLines <- grep(annotationColName[1], info, ignore.case=TRUE) - 1

		if (is.null(sep)) {
		    ## Find out the separator (sep) by taking the first two line of data, and comparing them.
		    ##  we assume it is either "," or "\t".
	    	titleLine <- info[nMetaDataLines + 1]
			sepNum <- regexpr('\t', titleLine)
			if (sepNum >= 2) {
				sep <- '\t'
			} else {
				sepNum <- regexpr(',', titleLine)
				if (sepNum >= 2) {
					sep <- ','
				} else {
					stop('The seperator is not Tab or comma!\n Please sepecify the seperator used in the file!')
				}
			}
		}

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

		## Read in annotation data
		annotation <- read.table(annotationFile, sep=sep, colClasses="character", header=TRUE, skip=nMetaDataLines,
		 	blank.lines.skip=TRUE, row.names=NULL, check.names=FALSE, quote=quote, comment.char="", strip.white=TRUE, fill=TRUE)

		colnames(annotation) <- toupper(colnames(annotation))
		## Create unique Id based on 50mer sequence
		nuID <- sapply(annotation[, toupper(annotationColName['sequence'])], seq2id)
		## check the TargetID first
		ann_target <- annotation[, toupper(annotationColName['target'])]
		comm_target <- id[id %in% ann_target]
		if (length(comm_target) == 0) {
			## check the ProbeID if id does not match the TargetID
			ann_target <- annotation[, toupper(annotationColName['probe'])]
			comm_target <- id[id %in% ann_target]
			if (length(comm_target) == 0) {
				width <- nchar(ann_target[1])
				id <- formatC(as.numeric(id), width=width, flag='0', format='d')
				comm_target <- id[id %in% ann_target]
				if (length(comm_target) == 0) stop('The annotation file does not match the data!')
			}
		} 
		if (length(comm_target) != length(id)) {
			warning('The annotation file does not match the data. Partial ids cannot be replaced!')
		}
		names(nuID) <- ann_target

		newId <- id
		newId[id %in% ann_target] <- nuID[comm_target]
	} else if (!is.null(lib)) {
		if (require(lib, character.only=TRUE)) {
			## check the TargetID first
			newId <- mget(id, get(paste(lib, 'TARGETID2NUID', sep=''), mode='environment'), ifnotfound=NA)
			newId <- unlist(newId)
			if (length(which(!is.na(newId))) == 0) {
				usingTargetID <- FALSE
				## check the ProbeID if id does not match the TargetID
				newId <- mget(id, get(paste(lib, 'PROBEID2NUID', sep=''), mode='environment'), ifnotfound=NA)
				if (length(which(!is.na(newId))) == 0) {
					width <- nchar(ls(envir=get(paste(lib, 'PROBEID2NUID', sep=''), mode='environment'))[1])
					id <- formatC(as.numeric(id), width=width, flag='0', format='d')
					newId <- mget(id, get(paste(lib, 'PROBEID2NUID', sep=''), mode='environment'), ifnotfound=NA)
					if (length(which(!is.na(newId))) == 0) {
						targetID <- pData(featureData(x.lumi))$TargetID
						newId <- mget(targetID, get(paste(lib, 'TARGETID2NUID', sep=''), mode='environment'), ifnotfound=NA)
						if (length(which(!is.na(newId))) == 0) stop('The library does not match the data!')
					}
				} 
			} else {
				usingTargetID <- TRUE
			}
			## Check for the targetIDs cannot be found in the lib.
			## Some known control genes will not be checked.
			naInd <- is.na(newId)
			controlId <- c('lysA','pheA','thrB','trpF', 'bla1','bla2','cat1','cat2','cre1','cre2','e1a1',
			'e1a2','gfp1','gfp2','gst1','gst2','gus1','gus2','lux1','lux2','lysA','neo1',
			'neo2','pheA','thrB','trpF')
			if (!usingTargetID) {
				TargetID <- featureData(x.lumi)$TargetID
				if (is.null(TargetID)) {
					if (!all(TargetID[naInd] %in% controlId)) {
						if (length(which(naInd)) < 10) {
							warning(paste('Identifiers:', paste(TargetID[naInd], collapse=','), ' cannot be found in the ', lib, '!', sep=''))
						} else {
							warning(paste('Some identifiers cannot be found in the ', lib, '!', sep=''))
						}
					}
				}
			} else if (!all(id[naInd] %in% controlId)) {
				if (length(which(naInd)) < 10) {
					warning(paste('Identifiers:', paste(id[naInd], collapse=','), ' cannot be found in the ', lib, '!', sep=''))
				} else {
					warning(paste('Some identifiers cannot be found in the ', lib, '!', sep=''))
				}
			}
			newId[naInd] <- id[naInd]
		} else {
			stop(paste('Annotation library', lib, 'is not installed!'))
		}
	} else {
		annotation <- pData(featureData(x.lumi))
		names(annotation) <- toupper(names(annotation))
		if (!is.null(annotation)) {
			sequence <- annotation[, 'PROBE_SEQUENCE']
			if (!is.null(sequence)) {
				cat('Directly converting probe sequence to nuIDs ...')
				newId <- sapply(sequence, seq2id)
				names(newId) <- id				
			} else {
				warning('Please provide the annotation file or lumi annotation library!')
			}
		} else {
			warning('Please provide the annotation file or lumi annotation library!')
		}
	}
	if (all(newId == id)) {
		conversion <- FALSE
	} else {
		conversion <- TRUE
	}

	if (any(duplicated(newId)))  {
		warning('Duplicated IDs found and were merged!')
		dupId <- unique(newId[duplicated(newId)])
		## determine whether the detection p-value close to 0 or 1 is significant
		detect.low <- exprs[which.max(detection(x.lumi)[,1]), 1]
		detect.high <- exprs[which.min(detection(x.lumi)[,1]), 1]
		
		rmIndex <- NULL
		for (dupId.i in dupId) {
			dupIndex <- which(newId == dupId.i)
			ave.exp <- colMeans(exprs(x.lumi)[dupIndex, ])
			exprs(x.lumi)[dupIndex[1],] <- ave.exp
			if (is(x.lumi, 'LumiBatch')) {
				totalBeadNum <- colSums(beadNum(x.lumi)[dupIndex, ])
				if (detect.low < detect.high) {
					maxDetection <- apply(detection(x.lumi), 2, min)
				} else {
					maxDetection <- apply(detection(x.lumi), 2, max)
				}

				temp <- colSums(se.exprs(x.lumi)[dupIndex,]^2 * (beadNum(x.lumi)[dupIndex,] - 1))
				temp <- temp / (totalBeadNum - length(dupIndex))
				se.exprs(x.lumi)[dupIndex[1],] <- sqrt(temp * (colSums(1/beadNum(x.lumi)[dupIndex,])))
				detection(x.lumi)[dupIndex[1],] <- maxDetection
				beadNum(x.lumi)[dupIndex[1],] <- totalBeadNum
			}
			rmIndex <- c(rmIndex, dupIndex[-1])
		}
		## remove duplicated
		x.lumi <- x.lumi[-rmIndex, ]
		newId <- newId[-rmIndex]
	}

	## update the feature names (probe ids)
	if (conversion) {
		names(newId) <- NULL
		featureNames(x.lumi) <- newId
		## update the feautre data
		featureData <- featureData(x.lumi)
		rownames(pData(featureData)) <- newId
		if (!is.null(pData(featureData)[,'PROBE_SEQUENCE'])) pData(featureData)[,'PROBE_SEQUENCE'] <- NULL
		featureData(x.lumi) <- featureData
		if (!is.null(lib)) annotation(x.lumi) <- lib

		## Add history tracking
		if (is(x.lumi, 'LumiBatch')) {
			history.finished <- as.character(Sys.time())
			history.command <- capture.output(print(match.call(addNuId2lumi)))
			if (is.null(x.lumi@history$lumiVersion)) x.lumi@history$lumiVersion <- rep(NA, nrow(x.lumi@history))
			lumiVersion <- packageDescription('lumi')$Version
			x.lumi@history<- rbind(x.lumi@history, data.frame(submitted=history.submitted, 
					finished=history.finished, command=history.command, lumiVersion=lumiVersion))
		}
	}

	return(x.lumi)
}