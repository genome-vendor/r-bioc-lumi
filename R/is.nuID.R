`is.nuID` <-
function(id) {

	if (length(id) > 1) {
		return(sapply(id, is.nuID))
	} else {
		if (nchar(id) < 2) return(FALSE)
		code <- c(LETTERS, letters, as.character(0:9), '_', '.')
		ind <- 1:length(code)
		names(ind) <- code
		if(length(grep('[^a-zA-Z0-9_.]', id)) > 0)  return(FALSE)

		id <- substring(id, 1:nchar(id), 1:nchar(id))
		num <- as.numeric(ind[id]) - 1
		checkCode <- num[1]
		num <- num[-1]
		if (checkCode == 63) return(FALSE)

		cutLen <- checkCode %% 3
		res <- floor(checkCode/3)

		codon <- rbind(floor(num / 4^2), floor((num %% 4^2) / 4), num %% 4)

		## Check whether the checkCode matches the sequence
		checkSum <- sum(num)
		if (res != checkSum %% 21) {
			return(FALSE)
		} else {
			return(TRUE)
		}
	}
}

