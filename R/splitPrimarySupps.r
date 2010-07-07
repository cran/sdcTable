splitPrimarySupps <- function(fullData) {
	upperInd <- function(index, dim) {
		ind <- which(!is.na(as.integer(do.call(rbind, lapply(dim, function(x) { match(index, x) } )))))[1]
		return(sort(dim[[ind]])[1])
	}	
	fdsub <- fullData$data[fullData$supps2check,]
	tmp <- list()
	for(i in 1:fullData$numberindexvars)
		tmp[[i]] <- sapply(fdsub[,fullData$indexvars[i]], upperInd, dim=fullData$dimensions[[i]])
	tmp <- as.data.frame(tmp)
	colnames(tmp) <- rep(" ", fullData$numberindexvars)
	
	f <- apply(tmp, 1, function(x) { paste(x, collapse="") } )
	fdsub <- cbind(fdsub, tmp)
	spl <- split(fdsub, f)
	spl	
}