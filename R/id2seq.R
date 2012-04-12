`id2seq` <-
function(id) {
	if (length(id) > 1) {
		return(sapply(id, id2seq))
	} else {
		## the coding is based on Base64 Content-Transfer-Encoding format (http://www.freesoft.org/CIE/RFC/1521/7.htm)
		#code <- c(LETTERS, letters, as.character(0:9), '+', '/')
		code <- c(LETTERS, letters, as.character(0:9), '_', '.')
		ind <- 1:length(code)
		names(ind) <- code
		if(length(grep('[^a-zA-Z0-9_.]', id)) > 0)  stop('Input id is not a nuID!')

		id <- substring(id, 1:nchar(id), 1:nchar(id))
		num <- as.numeric(ind[id]) - 1
		checkCode <- num[1]
		num <- num[-1]
		if (checkCode == 63)  warning('Coding error or not a nuID!\n Check code should not include "."!')
		cutLen <- checkCode %% 3
		res <- floor(checkCode/3)

		codon <- rbind(floor(num / 4^2), floor((num %% 4^2) / 4), num %% 4)

		## Check whether the checkCode matches the sequence
		checkSum <- sum(num)
		if (res != checkSum %% 21) {
			warning('Coding error or not a nuID!')
		}

		## nucleotide
		nucleotide <- c('A', 'C', 'G', 'T')
		codon <- nucleotide[codon + 1]
		len <- length(codon)
		seq <- paste(codon[1:(len-cutLen)], collapse='')
		return(seq)
	}
}

