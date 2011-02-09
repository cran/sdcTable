cellInfo <- function(outObj, characteristics, varNames, return=FALSE) {
	if ( !class(outObj) %in% c("outObj", "safeTable") )
		stop("the input object needs to be of class 'outObj' or 'safeTable'!\n")
	
	# calculate cell-strID from input (may be used in suppressCells() too!
	if ( class(outObj) == "safeTable" )
		nrDimVars <- length(outObj$outObj)-12		
	if ( class(outObj) == "outObj" )
		nrDimVars <- length(outObj$strObj$varNames)
	
	if ( length(characteristics) != nrDimVars | length(varNames) != nrDimVars )
		stop("please check your input objects!\n")			

	problematicPrimaryCells <- primSupp <- secondSupp <- NULL
	if ( class(outObj) == "outObj" ) {		
		str <- NULL
		for ( i in 1:nrDimVars ) {
			indexVar <- match(varNames[i], outObj$strObj$varNames)			
			str <- paste(str, outObj$levelObj[[indexVar]]$codesStandard[outObj$levelObj[[indexVar]]$codesOrig==characteristics[i]], sep="")
		}
		index <- match(str, outObj$fullTabObj$strID)
		if ( length(index) != 1 | is.na(index) )
			stop("Something has gone wrong. Please check your input parameter!\n")
		
		ID <- outObj$fullTabObj$strID[index]	
		cellStatus <- outObj$fullTabObj$status[index]

		if ( cellStatus == "u" ) {
			primSupp <- TRUE
			secondSupp <- FALSE
			if ( return == FALSE )
				cat ("The cell is a sensitive cell!\n")		
		}				
		if ( cellStatus == "s" ) {
			primSupp <- FALSE
			secondSupp <- FALSE				
			if ( return == FALSE )
				cat ("The cell is a candidate to be secondary suppressed!\n")
		}
		if ( cellStatus == "z" ) {
			primSupp <- FALSE
			secondSupp <- FALSE		
			if ( return == FALSE )
				cat ("The cell will be enforced for publication!\n")			
		}		
	}
	if ( class(outObj) == "safeTable" ) {
		returnObj <- NULL
		suppsInfo <- outObj$suppsInfo
		outObj <- outObj$outObj
		
		index <- 1:length(outObj$strID)
		for (i in 1:nrDimVars)
			index <- intersect(index, which(outObj[[varNames[i]]] == characteristics[i]))
		
		if (length(index) != 1)
			stop("something is completely wrong!\n")
		
		# output useful information	
		ID <- outObj$strID[index]	
		if ( outObj$status[index] %in% c("s", "z") ) {	
			primSupp <- FALSE
			secondSupp <- FALSE
			if ( return == FALSE )
				cat ("This cell can be published!\n")
		}
		if ( outObj$status[index] == "u" ) {
			primSupp <- TRUE		
			secondSupp <- FALSE		
			if ( return == FALSE ) {
				cat ("The cell with given characteristics cell has been marked as primary sensitive!\n")
			}
		}
		if ( outObj$status[index] == "x" ) {
			primSupp <- FALSE		
			secondSupp <- TRUE	
			sSupps <- sapply(suppsInfo, function(x) { x$secondSupps }) 
			pSupps <- sapply(suppsInfo, function(x) { x$primSupps }) 
			
			erg <- which(!is.na(sapply(sSupps, function(x) { match(outObj$strID[index], x)})))
			problematicPrimaryCells <- list()
			if ( length(erg) > 0 ) {
				res <- pSupps[c(erg)]	
				for ( i in 1:length(res) ) {
					out <- matrix(NA, nrow=length(res[[i]]), ncol=nrDimVars)
					for ( j in 1:length(res[[i]]) ) {
						for ( z in 1:nrDimVars ) {
							cellIndex <- match(res[[i]][j], outObj$strID)
							out[j,z] <- outObj[[5+z]][cellIndex]						
						}				
					}
					
					varNamesOrig <- names(outObj[13:length(outObj)])
					out <- out[,match(varNames, varNamesOrig)]			
					
					colnames(out) <- varNames
					out <- as.data.frame(out)
					problematicPrimaryCells[[i]] <- out
					if ( return == FALSE ) {
						cat ("The cell with given characteristics needs to be suppressed due to any of the following primary suppressed cells:\n")
						cat("SubTable", i,":\n")
						print(as.data.frame(out))					
					}
				}
			}
			else 
				stop("something went wrong!\n")
		}		
	}	
	if ( return == TRUE )
		return(list(ID=ID, characteristics=characteristics, varNames=varNames, primSupp=primSupp, secondSupp=secondSupp, problematicPrimaryCells=problematicPrimaryCells))
}