`seq2id` <-
function(seq) {

	if (length(seq) > 1) {
		return(sapply(seq, seq2id))
	} else {
		## the coding is based on Base64 Content-Transfer-Encoding format (http://www.freesoft.org/CIE/RFC/1521/7.htm)
		##code <- c(LETTERS, letters, as.character(0:9), '+', '/')
		code <- c(LETTERS, letters, as.character(0:9), '_', '.')
		seq <- toupper(seq)

		if(length(grep('[^ACGTU]', seq)) > 0)  stop('Input sequence is not a nucleotide sequnce!')

		## translate A, C, G, T, U, as 0, 1, 2, 3, 3
		num <- chartr("ACGTU", "01233", seq)

		## break the sequence as individual characters
		num <- substring(num, 1:nchar(num), 1:nchar(num))
		num <- as.numeric(num)
		if (length(num) %% 3 != 0) {
			appLen <- 3 - length(num) %% 3
			num <- c(num, rep(0, appLen))
		} else {
			appLen <- 0
		}

		codon <- matrix(num, nrow=3)
		code64 <- codon[1,] * 4^2 + codon[2,] * 4 + codon[3,]

		## check sum error checking
		checkSum <- sum(code64)
		res <- checkSum %% 21

		checkCode <- code[res*3 + appLen + 1]
		id <- paste(checkCode, paste(code[code64 + 1], collapse=''), sep='')
		return(id)
	}
}

