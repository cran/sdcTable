st <- function(subTab, method="orig") {
	xx <- as.data.frame(subTab)
	xx$Freq2 <- xx$Freq
	xx$Freq2[xx$status=="u"] <- -1
	xx$Freq2[xx$status=="x"] <- -2
	xx$Freq2[xx$status=="z"] <- 0
	
	xx$V1 <- substr(xx$strID, 1, 7)
	xx$V2 <- substr(xx$strID, 8, 13)
	
	if ( method=="orig")
		yy <- as.matrix(xtabs(xx$Freq~xx$V1+xx$V2))
	else
		yy <- as.matrix(xtabs(xx$Freq2~xx$V1+xx$V2))
	yy <- yy[-which(apply(yy,1,sum)==0),]
	yy
}


protectHYPERCUBE <- function(subTab, levelObj, strObj, allowZeros, protectionLevel, suppMethod, debug=TRUE) {
	# recode Index-Variables according to...? [DONE]
	f.recodeIndexVars <- function(subTabObj, levelObj, strObj) {
		#splitIDs <- sapply(subTabObj$strID, splitStrID, strInfo)
		splitIDs <- sapply(subTabObj$strID, splitStrID, strObj$strInfo)
		indexList <- list()
		for ( i in 1:length(levelObj) ) {
			a <- data.frame(a=splitIDs[i,])
			b <- data.frame(a=sort(as.character(unique(a$a))), c=1:length(unique(a$a)))
			indexList[[i]] <- as.character(b$a[match(a$a, b$a)])
			#indexList[[i]] <- b$c[match(a$a, b$a)]-1
		}
		indexList
	}
	
	# for a given supressed value, calculate a list of all 
	# possible diametral indices [DONE]
	f.diametralIndices <- function(subTabObj, indexList, cellToProtect) {		
		diametralIndices <- do.call("cbind", indexList)
		indToProtect <- diametralIndices[cellToProtect, ]
		
		for ( i in 1:ncol(diametralIndices) ) {
			if ( length(unique(diametralIndices[,i])) > 1 )
				diametralIndices <- subset(diametralIndices, diametralIndices[,i] != indToProtect[i])		
		}
		
		# is the status of any diametral indices == "z"? if so, remove!
		tmpDia <- apply(diametralIndices, 1, paste, collapse="")
		tmpIndex <- match(tmpDia, subTabObj$strID)
		zInd <- which(subTabObj$status[tmpIndex]=="z")
		zeroInd <- which(subTabObj$Freq[tmpIndex] == 0)
		
		if ( nrow(diametralIndices) == 0 )
			diametralIndices <- NULL
		# is it a primary or secondary suppression?
		SorP <- "P"
		#if ( subTabObj$status[cellToProtect] == "u" )
		#	SorP <- "P"
		if ( subTabObj$status[cellToProtect] == "x" )
			SorP <- "S"
		
		supp <- list()
		supp$cellToProtect <- cellToProtect
		supp$indToProtect <- indToProtect
		supp$SorP <- SorP
		supp$diametralIndices <- diametralIndices
		supp$zInd <- zInd
		supp$zeroInd <- zeroInd
		return(supp)
	}	
	# f.sub.calcInformationForSuppValg() calculates all the information needed for a quader which is given by g and its diametral value
	## subtab: a subtable of fullData$data containing suppressed values
	## supp: output object of diametralIndex(). containes the index of the actual value which needs to be protected and all possible diametral indices.
	## allowZeros, protectionLevel as above
	f.calcInformation <- function(subTabObj, diaObj, allowZeros, protectionLevel) {
		# g: cell to protect
		# d: a diametral cell
		calcQInfo <- function(g, d) {
			numberIndexVars <- length(g)
			
			# 1) identify quader
			l <- list()
			for ( i in 1:numberIndexVars )
				l[[i]] <- c(g[i], d[i])
			quader <- expand.grid(l)		
			quader <- quader[!duplicated(quader),]
			
			# 2) normquader
			dimQ <- nrow(quader)
			normQuader <- matrix(0, nrow=nrow(quader), ncol=ncol(quader))
			for ( i in 1:ncol(quader) ) 
				normQuader[which(quader[,i] != quader[dimQ,i]),i] <- 1
			normQuader		
			
			# 3) g|u indication?
			indexing <- rep("g", dimQ)
			indexing[which(apply(normQuader, 1, sum) %%2 != 0)] <- "u"
			return(list(quader=quader, normQuader=normQuader, indexing=indexing))
		}
		
		#nrDimensions <- length(diaObj$indToProtect)		
		# we create list which will contain the results for each diametral index of supp$indGeh
		resultObj <- list()
		
		# TODO: What is with 1-dimensional data?
		limit <- nrow(diaObj$diametralIndices)
		
		for ( z in 1:limit ) {
			# 1) we identify the quader given by supp$indGeh and supp$diametralIndices[z,]
			# TODO: What is with 1-dimensional data?
			qInfo <- calcQInfo(diaObj$indToProtect, diaObj$diametralIndices[z,])
			indexQuader <- qInfo$quader
			indexQuaderVec <- as.character(unlist(indexQuader))
			normQuader <- qInfo$normQuader
			indexing <- qInfo$indexing		
			
			# 2) position (indices==qPosition) of current quader in subTabObj
			valsQ <- pasteStrVec(as.character(unlist(indexQuader[,, drop=FALSE])), ncol(indexQuader))
			qPosition <- match(valsQ, subTabObj$strID)
			suppStatus <- subTabObj$status[qPosition]

			# 3) calculate various information about the selected quader (infoQuader)
			# 3.1) how many values would need to be suppressed for this quader
			indNonSupp <- which(subTabObj$status[qPosition] == "s" & subTabObj$status[qPosition] != "x")
			nrAdditionalSupps <- length(indNonSupp)
			
			# 3.2) whats the amount of information which needs to be suppressed?
			sumAdditionalSuppsFreq 	 <- sum(subTabObj$Freq[qPosition[indNonSupp]])
			sumAdditionalSuppsnumVal <- sum(subTabObj$numVal[qPosition[indNonSupp]])
			
			# 3.3) does the quader contains other single cells except for 
			# the primary suppressed value (diaObj$cellToPretect) to check?
			# subIndices = current quader without primary suppressed cell to check
			indSingleItems <- setdiff(which(subTabObj$Freq[qPosition]==1),1)
			singleItems <- NULL
			indikatorSingleItems <- FALSE
			if( length(indSingleItems) >= 1 ) {
				indikatorSingleItems <- TRUE
				singleItems <- indSingleItems
			}
			
			# 3.4) does the quader contain empty cells?
			indZeroItems <- setdiff(which(subTabObj$Freq[qPosition]==0),1)
			zeroItems <- NULL
			indikatorZeroItems <- FALSE		
			if( length(indZeroItems) > 0 ) {
				indikatorZeroItems <- TRUE
				zeroItems <- indZeroItems
			}	
			
			# 3.4.1) does the quader contain a cell that needs to be published?
			indProtected <- FALSE
			if ( any(subTabObj$status[setdiff(qPosition, qPosition[1])] == "z") )
				indProtected <- TRUE
			
			# 3.5) is the quader protected enough? (protectionLevel)
			# we need to check for interval-protection only if protectionLevel > 0
			schutzInd <- TRUE
			schutz <- protectionLevel
			if( protectionLevel > 0 ) {
				if( diaObj$SorP[1] == "P" ) {
					if ( !all(indexing=="u") )	{
						range <- min(subTabObj$Freq[qPosition[which(indexing=="u")]]) + min(subTabObj$Freq[qPosition[which(indexing=="g")]])
						X <- subTabObj$Freq[diaObj$cellToProtect]
						if( X == 0 ) {
							tmpInd <- which(subTabObj$status[qPosition] != "u" & subTabObj$Freq[qPosition] != 0)
							
							if( length(tmpInd) > 0 ) {
								# TODO: this needs testing !!! (page 60, repsilber)
								if( range <= min(subTabObj$Freq[tmpInd]) ) {
									schutzInd <- FALSE
									protectionLevel <- 0
								}
							}
						}
						else {
							schutz <- (100*range) / X
							if ( schutz < protectionLevel )
								schutzInd <- FALSE
						}
					}
				}
				
				
				# 4) return results
				# in this case, the cell is already protected, so we can stop!
				if( nrAdditionalSupps == 0 & schutzInd == TRUE & indikatorSingleItems == FALSE & indikatorZeroItems == FALSE ) {
					return(erg = NULL)
					break
				}
				
				resultObj[[z]] <- list(
						quaderStrID = valsQ,
						indexing = indexing,
						qPosition = qPosition,
						nrAdditionalSupps=nrAdditionalSupps,
						sumAdditionalSuppsFreq = sumAdditionalSuppsFreq,
						sumAdditionalSuppsnumVal = sumAdditionalSuppsnumVal,
						indikatorSingleItems = indikatorSingleItems,
						singleItems = singleItems,
						indikatorZeroItems = indikatorZeroItems,
						zeroItems = zeroItems,
						indikatorProtected = indProtected,
						schutz = schutz,
						schutzInd = schutzInd
				)				
			}
		}
		return(resultObj)
	}	
	
	# suppress a given quader
	suppressQuader <- function(subTabObj, qSupp) {		
		suppIndex <- setdiff(qSupp$qPosition, qSupp$qPosition[which(subTabObj$status[qSupp$qPosition]=="u")])
		subTabObj$status[suppIndex] <- "x"
		return(subTabObj)
	}	
	
	f.selectQuader <- function(subTabObj, diaObject, cellToProtect, allowZeros, suppMethod, protectionLevel, debug) {
		if ( !suppMethod %in% c("minSupps", "minSum", "minSumLogs" ) )
			stop("Error: not a valid suppression method!\n")
		
		if ( allowZeros==FALSE & sum(diaObject$zInd, diaObject$zeroInd) > 0 ) 			
			diaObject$diametralIndices <- diaObject$diametralIndices[-unique(c(diaObject$zInd, diaObject$zeroInd)),,drop=FALSE]
		
		infoObj <- f.calcInformation(subTabObj, diaObject, allowZeros, protectionLevel)
		
		# already protected
		if ( is.null(infoObj) ) {
			suppObj <- NULL
		}
		else {
			# which elements of iqsInfo are NULL?
			nullElements <- which(unlist(lapply(lapply(infoObj, '[[', 'qPosition'), function(x) { length(x) } )) == 0)
			if ( length(nullElements) > 0 )
				infoObj <- infoObj[-nullElements]
			
			# remove quaders that contain cells that need to be published
			protectedElements <- as.logical(do.call(rbind, lapply(infoObj, '[', 'indikatorProtected')))
			protectedElements <- which(protectedElements == TRUE)
			if ( length(protectedElements) > 0 )
				infoObj <- infoObj[-protectedElements]	
			
			if ( length(infoObj) == 0 ) {
				suppObj <- NA
				return(suppObj)
			}				
			# put iqs together so that we can choose the optimal suppression scheme
			qIndexNr <- 1:length(infoObj)
			nrAdditionalSupps <- as.numeric(as.character(do.call(rbind, lapply(infoObj, '[', 'nrAdditionalSupps'))))
			sumAdditionalSuppsFreq <- as.numeric(as.character(do.call(rbind, lapply(infoObj, '[', 'sumAdditionalSuppsFreq'))))
			indikatorSingleItems <- as.logical(do.call(rbind, lapply(infoObj, '[', 'indikatorSingleItems')))
			indikatorZeroItems <- as.logical(do.call(rbind, lapply(infoObj, '[', 'indikatorZeroItems')))
			schutz <- as.numeric(as.character(do.call(rbind, lapply(infoObj, '[', 'schutz'))))
			schutzInd <- as.logical(do.call(rbind, lapply(infoObj, '[', 'schutzInd')))
			
			possQuaders <- data.frame(qIndexNr, nrAdditionalSupps, sumAdditionalSuppsFreq, indikatorSingleItems, indikatorZeroItems, schutz, schutzInd)
			
			# all rows are removed where indicatorZeroItems != FALSE if allowZeros is FALSE
			if( allowZeros == FALSE ) {
				if ( !all(possQuaders$indikatorZeroItems == TRUE) )
					possQuaders <- possQuaders[possQuaders$indikatorZeroItems == FALSE,,drop=FALSE]	
				else {
					print(str(diaObj))
					stop("No pattern with empty-cells available!\n")
				}					
			}
			
			# are there any suppression schemes satisfying the necessary interval protection?
			indexIntervallOk <- FALSE
			if ( any(possQuaders$schutzInd==TRUE) ) 
				indexIntervallOk <- TRUE
			
			# do suppression schemes exist that do not contain single values? 
			# these are preferred suppression schemes.
			existNonSingles <- FALSE
			if ( any(possQuaders$indikatorSingleItems == FALSE) ) {
				existNonSingles <- TRUE		
				possQuaders <- possQuaders[possQuaders$indikatorSingleItems==FALSE,,drop=FALSE]
			}	
			
			if( indexIntervallOk == TRUE ) {
				if( min(possQuaders$nrAdditionalSupps) > 0 & debug == TRUE)
					cat("# additional secondary Supps:", min(possQuaders$nrAdditionalSupps)," ")
				if( suppMethod == "minSupps" ) {
					possQuaders <- possQuaders[which(possQuaders$nrAdditionalSupps == min(possQuaders$nrAdditionalSupps)),,drop=FALSE]
					possQuaders <- possQuaders[which(possQuaders$sumAdditionalSuppsFreq == min(possQuaders$sumAdditionalSuppsFreq)),,drop=FALSE]
				}
				if( suppMethod == "minSum" ) {
					possQuaders <- possQuaders[which(possQuaders$sumAdditionalSuppsFreq == min(possQuaders$sumAdditionalSuppsFreq)),,drop=FALSE]
					possQuaders <- possQuaders[which(possQuaders$nrAdditionalSupps == min(possQuaders$nrAdditionalSupps)),,drop=FALSE]
				}
				if( suppMethod == "minSumLogs" ) {
					possQuaders <- possQuaders[which(possQuaders$sumAdditionalSuppsFreq == min(log(1+possQuaders$sumAdditionalSuppsFreq))),,drop=FALSE]
					possQuaders <- possQuaders[which(possQuaders$nrAdditionalSupps == min(possQuaders$nrAdditionalSupps)),,drop=FALSE]
				}
				
				# finally choose the suppression scheme
				possQuaders <- possQuaders[1,]		
				suppObj <- infoObj[[possQuaders$qIndexNr]]		
			}
			# we have a problem: 
			# there is no suppression scheme satisfying the nessessary interval protection
			else {
				# all cells in this subtable are either primary or secondary suppressed
				# in this case everything is ok and we have no problem!
				if( all(subTabObj$status=="u" | subTabObj$status == "x") ) 
					suppObj <- NULL
				# no suppression scheme satisfies the required interval protection
				# the suppression pattern with the max. protection level is selected 
				else {
					possQuaders <- possQuaders[which(possQuaders$schutz == max(possQuaders$schutz)),]
					possQuaders <- possQuaders[which(possQuaders$nrAdditionalSupps == min(possQuaders$nrAdditionalSupps)),]
					possQuaders <- possQuaders[which(possQuaders$sumAdditionalSuppsFreq == min(possQuaders$sumAdditionalSuppsFreq)),]
					possQuaders <- possQuaders[1,]
					suppObj <- infoObj[[possQuaders$qIndexNr]]
				}
			}			
		}
		suppObj
	}
			
	f.selectAdditionalQuader <- function(subTabObj, levelObj, strObj, suppObject, diaObject) {
		fn.indices <- function(g,d) {
			l <- list()
			numberIndexVars <- length(g)
			for ( i in 1:numberIndexVars )
				l[[i]] <- c(g[i], d[i])
			quader <- expand.grid(l)		
			quader <- quader[!duplicated(quader),,drop=F]
			strIDs <- as.character(apply(quader, 1, paste, collapse=""))
			strIDs
		}		
		# is cellToProtect a marginal cell? if so -> no additional quader needs to be found		
		cellInfo <- isMarginalSum(subTabObj$strID, strObj$strInfo) 
		
		# calculate complete indices
		tmpStr <- apply(diaObject$diametralIndices, 1, function(x) { fn.indices(diaObject$indToProtect,x)})
		indices <- apply(tmpStr, 1, function(x) { match(x, subTabObj$strID) } )
		
		# remove suppObject from diaObject		
		if ( length(suppObject)==1 )
			suppObject <- suppObject[[1]]
		indexSuppObj <- which(apply(indices, 1, function(x) { all(suppObject$qPosition %in% x)}))
		if ( length(indexSuppObj)!=1 )
			stop("Error!\n")
		
		indices <- indices[-indexSuppObj,]
		diaObject$diametralIndices <- diaObject$diametralIndices[-indexSuppObj,,drop=FALSE]
		
		cellToProtect <- diaObject$cellToProtect
		warn <- FALSE
		if ( !(cellToProtect %in% cellInfo$indexTotCells) ) {	
			# if we are dealing with a table that includes only one 
			# respondent (all marginals == one inner cell) we do not
			# have to find another quader
			# 0 0 0 | 0
			# 0 1 0 | 1
			# 0 0 0 | 0
			# ------|
			# 0 1 0 | 1					
			if ( length(which(subTabObj$Freq[cellInfo$indexInnerCells]!=0 )) > 1 ) {
				indexSingles <- suppObject$qPosition[suppObject$singleItems]
				
				# if the remaining singleton is in a margin,	
				# it is not required to find an additional quader
				# if none of the indexSingles are marginal cells, then we have to find additional suppressions
				if ( !any(indexSingles %in% cellInfo$indexTotCells) ) {
					# task: find a quader that does not contain indexSingles as qPosition
					tt <- diaObject
					workIndex <- indices
					indexRem1 <- which(apply(indices, 1, function(x) { any(indexSingles %in% x )}))
					indexRem2 <- which(apply(indices, 1, function(x) { any(subTabObj$status[x]=="z")}))
					indexRem3 <- which(apply(indices, 1, function(x) { any(subTabObj$Freq[x]==0)}))
					indexRem <- sort(unique(c(indexRem1, indexRem2, indexRem3)))
					
					if ( length(indexRem) == nrow(indices) )  {
						suppObj <- NULL
						return(suppObj)
					}											
					
					workIndex <- workIndex[-indexRem,,drop=F]
					tt$diametralIndices <- tt$diametralIndices[-indexRem,,drop=FALSE]
					infoObj <- f.calcInformation(subTabObj, tt, allowZeros, protectionLevel)

					# find quader with least additional supps (-> newly calculated)
					if ( length(infoObj) > 0 ) {
						newSupps <- rep(NA, length(infoObj))
						totNr <- length(infoObj[[1]]$qPosition)
						for (i in 1:length(infoObj))
							newSupps[i] <- totNr - length(which(subTabObj$status[infoObj[[i]]$qPosition] %in% c("u","x")))
						infoObj <- infoObj[which(newSupps==min(newSupps))]
						if ( length(infoObj) == 1 ) 
							suppObj <- infoObj
						else {
							suppObj <- infoObj[which.min(lapply(infoObj, function(x) { x$nrAdditionalSupps } ))]
						}
					}
					# no additional cube could be found
					else {
						suppObj <- NULL
					}						
				}
				# singletons in margins -> no additional quader required
				else {
					suppObj <- suppObject
				}
			}
			# only one cell != 0 in the inner cells of the table	
			else {				
				suppObj <- suppObject
			}
		}		
		return(suppObj)
	}		
	
	######################
	### start programm ###
	######################
	
	# prepare subTabObj by recoding index-variables
	indexList <- f.recodeIndexVars(subTab, levelObj, strObj)
	cellInfo <- isMarginalSum(subTab$strID, strObj$strInfo) 
	cellsToProtect <- which(subTab$status %in% c("u"))
	nrSupps <- nrSuppStart <- length(cellsToProtect)
	
	# protect all primary suppressed cells sequentially
	if ( debug == TRUE )
		cat("\nThe iterative protection procedure is now started ... \n")
	#for ( i in 1:7) {
	for ( i in 1:length(cellsToProtect) ) {
		if ( debug == TRUE )
			cat("--> Cell",i,"|",length(cellsToProtect)," (ID:",subTab$strID[cellsToProtect[i]],")...")
		
		# calculate diametral indices for current cell g=cellsToProtect[i]
		# -> diaObj: all diametral indices and information about the suppressed value 
		diaObj <- f.diametralIndices(subTab, indexList, cellsToProtect[i]) 
		
		# if no cube can be found that does not contain cells with status "z"
		# we stop and have to relax!	
		if ( length(diaObj$zInd) == nrow(diaObj$diametralIndices) ) {
			subTab <- NULL
			break
		}
		# we look for an "optimal" cube and suppress additional cells 
		# if it is necessary
		suppObj <- f.selectQuader(subTab, diaObj, cellsToProtect[i], allowZeros, suppMethod, protectionLevel, debug)
		if ( length(suppObj)!=0 && is.na(suppObj) ) {
			subTab <- NULL
			break			
		}
		if ( !is.null(suppObj) && !is.na(suppObj) ) {
			if ( any (subTab$Freq[suppObj$qPosition] == 0) )
				stop("Cells with Frequency=0 had to be suppressed in selectQuader()\n")			
			subTab <- suppressQuader(subTab, suppObj)
			
			# additional quader needs to be found 
			# only if it is not a single value in the margins
			# and the cube includes cells with frequency=1
			indToContinue <- !(cellsToProtect[i] %in% cellInfo$indexTotCells & subTab$Freq[cellsToProtect[i]] == 1)
			if ( suppObj$indikatorSingleItems==TRUE & indToContinue==TRUE ) {
				# find additional cube that does not contain the single cells
				suppObj <- f.selectAdditionalQuader(subTab, levelObj, strObj, suppObj, diaObj)
				if ( !is.null(suppObj) ) {
					subTab <- suppressQuader(subTab, suppObj)
					if ( any(subTab$Freq[suppObj$qPosition] == 0) )
						stop("Cells with Frequency=0 had to be suppressed in selectAdditionalQuader()\n")
				}					
				else {
					#cat("wir haben ein Problem!\n")
					subTab <- NULL
					break
				}								
			}				
		}
		if ( debug == TRUE )
			cat("[DONE]\n")			
	}	
	return(subTab)
}
