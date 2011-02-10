# split the complete string of level-information according to  levelObj
splitStrVec <- function(strObj, levelObj) {
	vars <- match(unlist(lapply(levelObj, function(x) { x$varName } )), strObj$varNames)
	v <- strObj$strID
	info <- strObj$strInfo[vars]	
	if( length(vars) == 1 ) 
		newStr <- as.vector(sapply(v, substr, info[[1]][1], info[[1]][2]))		
	else { 
		colInd <- unlist(lapply(info, function(x) { seq(x[1], x[2])}))	
		newStr <- matrix(unlist(strsplit(v,"")), nrow=length(v), byrow=T)[,colInd]
		newStr <- apply(newStr, 1, paste, collapse="")		
	}	
	strObj <- list()
	strObj$strID <- newStr
	strObj$strInfo <- info
	strObj	
}

## Wrapper function for pasting key-Variables
pasteStrVec <- function(strVec, nrVars) {
	if(length(strVec) %% nrVars != 0)
		stop("Wrong Dimensions!\n")
	else {
		.Call( "myPaste", as.character(strVec), nrVars, PACKAGE = "sdcTable")
	}
}

# alternative to expand.grid (used for pasteStrVec!)
expand <- function(inputList, vector=TRUE) {
	uniques <- sapply(inputList, length)
	nrPoss <- prod(uniques)
	if ( vector == TRUE ) {
		out <- NULL
		for ( i in 1:length(inputList) ) {
			if ( i == 1 ) 
				out <- rep(inputList[[i]], nrPoss/length(inputList[[i]]))
			else 
				out <- c(out, rep(inputList[[i]], each=prod(uniques[1:(i-1)]), nrPoss/length(rep(inputList[[i]], each=prod(uniques[1:(i-1)])))))
		}		
	}
	else {
		out <- list()
		for ( i in 1:length(inputList) ) {
			if ( i == 1 ) 
				out[[i]] <- rep(inputList[[i]], nrPoss/length(inputList[[i]]))
			else 
				out[[i]] <- rep(inputList[[i]], each=prod(uniques[1:(i-1)]), nrPoss/length(rep(inputList[[i]], each=prod(uniques[1:(i-1)]))))
		}		
	}
	out
}	

# split stringID's according to it's strInfo
splitStrID <- function(strId, strInfo) {
	out <- c()
	for ( i in 1:length(strInfo) )
		out <- c(out, substr(strId, strInfo[[i]][1], strInfo[[i]][2]))
	out
}	

# the same as splitStrID() with the difference being the output format (list)
splitStrID2 <- function(strId, strInfo) {
	out <- list()
	for ( i in 1:length(strInfo) )
		out[[i]] <- sort(unique(substr(strId, strInfo[[i]][1], strInfo[[i]][2])))
	out
}	

# returns all possible strIDs from a given level and dimObj
calcPossiblestrIDS <- function(currentLevel, dimObj) {
	dd <- list()
	for ( i in 1:length(dimObj) )
		dd[[i]] <- dimObj[[i]][[currentLevel[i]]]
	return(pasteStrVec(as.character(expand(dd)), length(dd)))	
} 	

# returns the strIDs and corresponding indices of inner and marginal table cells
isMarginalSum <- function(strIDs, strInfo) {		
	out <- splitStrID2(strIDs, strInfo)	
	innerCells <- apply(expand.grid(sapply(out, function(x) { x[2:length(x)] })), 1, paste, collapse="")
	totCells <- setdiff(strIDs, innerCells)
	indexTotCells <- match(totCells, strIDs)
	indexInnerCells <- match(innerCells, strIDs)
	return(list(innerCells=innerCells, totCells=totCells, indexInnerCells=indexInnerCells, indexTotCells=indexTotCells))
}
