changeCellStatus <- function(outObj, varNames, characteristics, rule, codesOrig=TRUE, suppZero=FALSE) {
	if (length(varNames) != length(characteristics) )
		stop("please check your input-objects!\n")
	
	if ( !rule %in% c("u","z","x","s") )
		stop("Please specify a correct rule!\n")
	
	# calculate strID of cell to suppress
	varsInLevelObj <- unlist(lapply(outObj$levelObj, function(x) { x$varName }))
	positionInLevelObj <- unlist(lapply(outObj$levelObj, function(x) { x$posIndex }))
	
	pos <- NULL
	strID <- NULL
	for ( i in 1:length(varNames) )	{
		pos[i] <- match(varNames[i], varsInLevelObj)
		if ( codesOrig == TRUE )
			strInd <- match(characteristics[i], outObj$levelObj[[pos[i]]]$codesOrig)
		else
			strInd <- match(characteristics[i], outObj$levelObj[[pos[i]]]$codesStandard)	
		strID[i] <- outObj$levelObj[[pos[i]]]$codesStandard[strInd]
	}	
		
	# variables available (in strID, levelObj,..) but not listed in varNames
	notUsedVars <- setdiff(varsInLevelObj, varNames)
	if ( length(notUsedVars) > 0 ) {
		for ( i in 1:length(notUsedVars) ) {
			pos <- c(pos, match(notUsedVars[i], varsInLevelObj))
			strID <- c(strID, outObj$levelObj[[pos[length(pos)]]]$codesStandard[1])
		}
	}	
	
	strID <- paste(strID[pos], collapse="")
	# suppress
	suppInd <- match(strID, outObj$fullTabObj$strID)
	if ( length(suppInd) != 1 )
		stop("something got wrong in suppressCells()\n")
	
	# set primary suppressed!
	if ( rule == "u" ) {
		if ( suppZero==FALSE & outObj$fullTabObj$Freq[suppInd] == 0 )
			cat("===> the given cell has been not been set as primary sensitive cell due to parameter 'suppZero'!\n")
		else {
			cat("===> the given cell has been set as primary sensitive cell!\n")
			outObj$fullTabObj$status[suppInd] <- "u"	
			outObj$fullTabObj$w[suppInd] <- 0
		}	
	}
	# mark cell as potential secondary suppression
	if ( rule == "s" ) {
		cat("===> the given cell will be marked as a candidate for secondary suppression. The previous state was: ''",outObj$fullTabObj$status[suppInd],"''.\n")
		outObj$fullTabObj$status[suppInd] <- "s"		
	}	
	# mark cell as secondary suppressed
	if ( rule == "x" ) {
		cat("===> the given cell will be marked as secondary suppressed. The previous state was: ''",outObj$fullTabObj$status[suppInd],"''.\n")
		outObj$fullTabObj$status[suppInd] <- "x"	
		outObj$fullTabObj$w[suppInd] <- outObj$fullTabObj$Freq[suppInd]
	}	
	# force cell to be published
	if ( rule == "z" ) {
		cat("===> the given cell is marked to be published! The previous state was: ''",outObj$fullTabObj$status[suppInd],"''.\n")
		outObj$fullTabObj$status[suppInd] <- "z"	
		outObj$fullTabObj$w[suppInd] <- outObj$fullTabObj$Freq[suppInd]
	}
	outObj
}
