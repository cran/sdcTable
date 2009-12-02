## calcSubset is a function to draw a specific subtable from the entire data-structure
## fullData: the entire data structure as given in Quaderverfahren() and qf()
## group: list-element from output of splitPrimarySupps()
calcSubset <- function(fullData, group) {
	calcIndices <- function(index, dim) {
		ind <- which(!is.na(as.integer(do.call(rbind, lapply(dim, function(x) { match(index, x) } )))))[1]
		return(dim[[ind]])
	}	
	tmp <- group[1,]
	indices <- list()
	for(i in 1:fullData$numberindexvars)
		indices[[i]] <- calcIndices(tmp[1,i], fullData$dimensions[[i]])
	
	subtab <- fullData
	for(i in 1:fullData$numberindexvars)
		subtab$data <- subset(subtab$data, subtab$data[,i] %in% indices[[i]])
	return(subtab)
}