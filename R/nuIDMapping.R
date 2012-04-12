nuID2probeID <- function(nuID, lib="lumiHumanV1") {
	if (length(nuID) == 0) return(NULL)

	if (!is.null(lib)) {
		if (length(grep('\\.db', lib)) > 0) {
			warning(paste(lib, 'does not include nuID conversion information!'))
			return(nuID)
		}
	}
	if (require(lib, character.only=TRUE)) {
		env <- get(paste(lib, 'PROBEID2NUID', sep = ""), mode = "environment")
		probe2nuID <- unlist(as.list(env))	
		allProbe <- names(probe2nuID)
		keepInd <- !is.na(probe2nuID)
		allProbe <- allProbe[keepInd]
		names(probe2nuID) <- allProbe
		probe <- lapply(nuID, function(x) allProbe[probe2nuID == x])
		names(probe) <- nuID
		return(probe)
	} else {
		print(paste(lib, ' annotation library is required!', sep=''))
	}
}

nuID2targetID <- function(nuID, lib="lumiHumanV1") {
	if (length(nuID) == 0) return(NULL)
	if (!is.null(lib)) {
		if (length(grep('\\.db', lib)) > 0) {
			warning(paste(lib, 'does not include nuID conversion information!'))
			return(nuID)
		}
	}
	if (require(lib, character.only=TRUE)) {
		env <- get(paste(lib, 'TARGETID2NUID', sep = ""), mode = "environment")
		target2nuID <- unlist(as.list(env))	
		allTarget <- names(target2nuID)
		keepInd <- !is.na(target2nuID)
		allTarget <- allTarget[keepInd]
		names(target2nuID) <- allTarget
		target <- lapply(nuID, function(x) allTarget[target2nuID == x])
		names(target) <- nuID
		return(target)
	} else {
		print(paste(lib, ' annotation library is required!', sep=''))
	}
}

probeID2nuID <- function(probeID, lib="lumiHumanV1") {
	if (length(probeID) == 0) return(NULL)
	if (!is.null(lib)) {
		if (length(grep('\\.db', lib)) > 0) {
			warning(paste(lib, 'does not include nuID conversion information!'))
			return(nuID)
		}
	}
	if (!require(annotate)) print('Please install "annotate" library!')
	if (require(lib, character.only=TRUE)) {
		nuID <- unlist(lookUp(probeID, lib, 'PROBEID2NUID'))
		return(nuID)
	} else {
		print(paste(lib, ' annotation library is required!', sep=''))
	}
}

targetID2nuID <- function(targetID, lib="lumiHumanV1") {
	if (length(targetID) == 0) return(NULL)
	if (!is.null(lib)) {
		if (length(grep('\\.db', lib)) > 0) {
			warning(paste(lib, 'does not include nuID conversion information!'))
			return(nuID)
		}
	}
	if (!require(annotate)) print('Please install "annotate" library!')
	if (require(lib, character.only=TRUE)) {
		nuID <- unlist(lookUp(targetID, lib, 'TARGETID2NUID'))
		return(nuID)
	} else {
		print(paste(lib, ' annotation library is required!', sep=''))
	}
}
