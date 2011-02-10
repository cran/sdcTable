#############################################
### calc a fullTabObj from given startObj ###
#############################################
calcFullTable <- function(dataset, levelObj, freqVar=NULL, numVar=NULL, weightVar=NULL, sampWeight=NULL ) {
	# prepare data prior to trying to anonymize input
	f.prepareInputData <- function(dat, levelObj, freqVar, numVar, weightVar, sampWeight) {	
		freqVarInd <- numVarInd <- weightVarInd <- sampWeightInd <- NA		
		if ( !is.null(freqVar) )
			freqVarInd <- match(freqVar, colnames(dat))	
		if ( !is.null(numVar) )
			numVarInd <- match(numVar, colnames(dat))
		if ( !is.null(weightVar) )
			weightVarInd <- match(weightVar, colnames(dat))			
		if ( !is.null(sampWeight) )
			sampWeightInd <- match(sampWeight, colnames(dat))				
		
		# quick fix for levelObj
		if ( !is.list(levelObj[[1]]) ) {
			lO <- list()
			lO[[1]] <- levelObj
			levelObj <- lO; rm(lO)
		}	
		
		posIndex <- sapply(levelObj, function(x) { x$posIndex } )
		
		# sort  levelObj according to the positions of the table variables
		# in the input-data
		levelObj <- levelObj[c(order(posIndex))]
		 
		varNames <- sapply(levelObj, function(x) { x$varName } )
		posIndex <- sapply(levelObj, function(x) { x$posIndex } )

		# check for NA's, ask to impute? remove?
		if ( any (is.na(dat[,posIndex])) ) {
			dimBefore <- dim(dat)
			NAs <- apply(dat[,posIndex], 2, function(x) { length(which(is.na(x))) } )
			cat("--> Info: there are missing values in your dimensional variables!\n")
			cat("--> Number of missings in each dimensional variable:\n")
			print(NAs)
			
			cat("--> Note: by default, we are ommiting these cases!\n")
			cat("--> Note: you could also impute them using R-package 'VIM'!\n")
			dat <- na.omit(dat)
			dimAfter <- dim(dat)
			cat("--> Note:", dimBefore[1] - dimAfter[1],"cases have been deleted!\n")
		}			
		
		if ( any(colnames(dat)[posIndex] != varNames) )
			stop ("please check input objects in prepareInputData()!\n")
		
		# remove or change levels if neccessary (levels with only 1 sublevel!)
		for ( i in 1:length(posIndex) ) {
			ll 	<- levelObj[c(i)][[1]]	
			# factors --> character
			dat[,ll$varName] <- as.character(dat[,ll$varName])
			
			# Fix: dat needs to have levelsOrig (not standard codes!)
			lDown <- ll$levelsRemoveOrig
			lUp	<- ll$levelsRemoveOrigUp
			if ( length(lDown) > 0 ) {
				# order by standard codes decreasingly!
				orderInd <- order(ll$codeRemoveOrig, decreasing=TRUE)
				lDown <- lDown[orderInd]
				lUp <- lUp[orderInd]					
				for ( j in 1:length(lDown) ) {
					index <- which(dat[,ll$varName] == lDown[j])
					# otherwise, the index has already been removed from the data
					if ( length(index) > 0 ) {
						if ( !is.null(freqVar) & !is.na(freqVarInd) )
							dat <- dat[-index,]
						else
							dat[index, ll$varName] <- lUp[j]
						#if ( !is.null(freqVar) & !is.na(freqVarInd) )
						#	dat[index, freqVarInd] <- rep(0, length(index))
					}							
				}
			}			
		}		
		
		# generate unique string for each cell | observation	
		strInfo <- list()
		for ( i in 1:length(levelObj) ) {
			sumCur <- sum(levelObj[[i]]$levelStructure)	
			if ( i == 1 )
				strInfo[[i]] <- c(1, sumCur)
			else 
				strInfo[[i]] <- c(1+max(strInfo[[c(i-1)]]), max(strInfo[[c(i-1)]])+sumCur)
		}
		
		# if non-standard input-codes are used, create strID with standardized codes
		tmpDat <- dat
		for ( i in 1:length(posIndex) ) {
			if ( all(unique(dat[,posIndex[i]]) %in% levelObj[[i]]$codesStandard) )
				tmpDat[,posIndex[i]] <- dat[,posIndex[i]]
			else
				tmpDat[,posIndex[i]] <- levelObj[[i]]$codesStandard[match(dat[,posIndex[i]], levelObj[[i]]$codesOrig)]
		}
			
		NAs <- apply(tmpDat[,posIndex,drop=FALSE], 2, function(x) { length(which(is.na(x))) } )
		colInd <- which(NAs > 0)
		if ( length(colInd) > 0 ) {
			for ( i in 1:length(colInd) )
				tmpDat <- tmpDat[-which(is.na(tmpDat[,colInd[i]])),]
		}
		dat <- tmpDat # remove if not working!
			
		strObj <- list()
		strObj$strID <- pasteStrVec(unlist(tmpDat[,posIndex, drop=FALSE]), length(posIndex))		
		strObj$strInfo <- strInfo
		strObj$varNames <- varNames
		
		datObj <- list()
		datObj$microData <- FALSE
		datObj$numVal <- rep(NA, length(strObj$strID))
		datObj$Freq <- rep(NA, length(strObj$strID))
		
		if ( !is.null(numVar) &  !is.na(numVarInd) )
			datObj$numVal <- dat[,numVarInd]				
		
		# microData!
		if ( is.null(freqVar) | is.na(freqVarInd) ) {
			datObj$microData <- TRUE
			datObj$Freq <- rep(1, length(strObj$strID))
			if ( !is.null(sampWeight) & !is.na(sampWeightInd) ) {
				datObj$Freq <- round(dat[,sampWeightInd])	
				if ( !is.null(numVar) & !is.na(numVarInd) )
					datObj$numVal <- round(dat[,numVarInd]*dat[,sampWeightInd])
			}				
		}		
		# data are already aggregated
		else if ( !is.null(freqVar) & !is.na(freqVarInd) ) {
			datObj$Freq <- dat[,freqVarInd]	
			if ( !is.null(sampWeight) & !is.na(sampWeightInd) ) {
				datObj$Freq <- round(datObj$Freq*dat[,sampWeightInd])	
				if ( !is.null(numVar) & !is.na(numVarInd) )
					datObj$numVal <- round(dat[,numVarInd]*dat[,sampWeightInd])
			}			
		}		
		else {
			stop("please check your input parameters!\n")
		}
		
		if ( is.null(weightVar) | is.na(weightVarInd) ) 
			datObj$w <- datObj$Freq
		else
			datObj$w <- dat[,weightVarInd]
		return(list(datObj=datObj, strObj=strObj, levelObj=levelObj))
	}	
		
	# generate a minimal table from microdata
	f.tableFromMicroData <- function(datObj, levelObj) {
		tableVars <- sapply(levelObj, function(x) x$varName)
		if( length(tableVars) != length(levelObj) )
			stop("please check your levelObj with respect to the names of the spanning variables in the dataset!\n")
		
		ag1 <- aggregate(datObj$Freq, by=list(strObj$strID), sum)
		ag2 <- aggregate(datObj$w, by=list(strObj$strID), sum)		
		ag3 <- aggregate(datObj$numVal, by=list(strObj$strID), sum)
		mintabObj <- list()
		mintabObj$strID <- ag1[,1]
		mintabObj$Freq <- ag1[,2]
		mintabObj$w <- ag2[,2]
		mintabObj$numVal <- ag3[,2]
		
		exDims <- sort(unique(strObj$strID)) 
		#possDims <- sort(pasteStrVec(as.character(expand(lapply(levelObj, function(x) { x$allDims[x$codesMinimal==TRUE]}), vector=TRUE)), length(levelObj)))
		possDims <- sort(pasteStrVec(as.character(expand(lapply(levelObj, function(x) { x$codesStandard[x$codesMinimal==TRUE]}), vector=TRUE)), length(levelObj)))
	
		# problematic are all levels that should exist, but do not exist
		# they are filled with NA/0
		probl <- setdiff(possDims, exDims)
		if ( length(probl) > 0 ) {
			mintabObj$strID <- c(mintabObj$strID, probl)
			mintabObj$Freq <- c(mintabObj$Freq, rep(0, length(probl)))
			mintabObj$w <- c(mintabObj$w, rep(0, length(probl)))
			mintabObj$numVal <- c(mintabObj$numVal, rep(0, length(probl)))
		}		
		mintabObj		
	}
	
	f.fillTable <- function(minTabObj, levelObj, strObj) {
		nrIndexvars <- length(levelObj)
		fullDims <- lapply(levelObj, function(x) {x$dimensions})
		
		# complete
		allDims <- pasteStrVec(expand(lapply(levelObj, function(x) x$codesStandard)), length(levelObj))
		
		# the subtotals that need to be calculated
		subTotals <- setdiff(allDims, minTabObj$strID)
		
		fullTabObj 			<- minTabObj
		fullTabObj$strID 	<- c(minTabObj$strID, subTotals)
		fullTabObj$Freq 	<- c(minTabObj$Freq, rep(NA, length(subTotals)))
		fullTabObj$w 		<- c(minTabObj$w, rep(NA, length(subTotals)))
		fullTabObj$numVal 	<- c(minTabObj$numVal, rep(NA, length(subTotals)))
		
		rI <- FALSE # run-Indicator
		# backup strObj <- we can edit strObjWork$strID later for performance reasons
		strObj$strID <- fullTabObj$strID
		strObj$indexVec <- 1:length(strObj$strID)
		strObjWork <- strObj
		
		while( rI == FALSE ) {
			for( i in 1:nrIndexvars ) {			
				cat("processing dimension",i,"|",nrIndexvars," ...\n"); flush.console()
				tt <- Sys.time()
				f1 <- splitStrVec(strObj, levelObj[c(i)])$strID
				
				for( j in (length(fullDims[[i]])):1 ) {					
					#cat("===> j=",j,"|",length(fullDims[[i]]),"\n")
					indRows <- f1 %in% fullDims[[i]][[j]]
					
					switch <- FALSE
					if ( all(indRows == TRUE) )
						switch <- TRUE
					
					rows <- fullTabObj$strID[indRows]
					strObjWork$strID <- strObj$strID[indRows]
					strObjWork$indexVec <- strObj$indexVec[indRows]
					# strID numeric for faster matching!
					strIDn <- as.numeric(strObjWork$strID)									
								
					if ( nrIndexvars > 1 ) {
						f <- splitStrVec(strObjWork, levelObj[c(-i)])$strID
						spl <- split(1:length(f), f)				
					}									
					else {
						spl <- list()
						spl[[1]] <- 1:length(rows)	
					}					
					
					for ( z in length(spl):1 ) {
						if ( switch == TRUE ) 
							posIndex <- spl[[z]]
						else	
							posIndex <- strObjWork$indexVec[match(as.numeric(rows[spl[[z]]]), strIDn)]

						freqs <- fullTabObj$Freq[posIndex]						
						ind <- which(is.na(as.numeric(freqs)))
						if( length(ind) == 1 ) {
							fullTabObj$Freq[posIndex[ind]] <- sum(as.integer(freqs), na.rm=T)
							nums <- fullTabObj$numVal[posIndex]
							ws <- fullTabObj$w[posIndex]
							fullTabObj$w[posIndex[ind]] <- sum(as.numeric(ws), na.rm=T)
							fullTabObj$numVal[posIndex[ind]] <- sum(as.numeric(nums), na.rm=T)
						}
					}	
				}
				cat("[DONE in:", round(Sys.time()-tt, 2),"]\n"); flush.console()
			}
			if( all(!is.na(fullTabObj$Freq)) ) 
				rI <- TRUE	
		}
		fullTabObj
	}	
	
	outObj <- list() # final output object
	
	cat("Time for f.prepareInputData():\n")
	print(system.time({						
	startObj <- f.prepareInputData(dataset, levelObj, freqVar, numVar, weightVar, sampWeight)			
	}))	
	
	datObj <- startObj$datObj
	strObj <- startObj$strObj
	levelObj <- startObj$levelObj
	
	# calc minimal table from microdata
	cat("Time for f.tableFromMicroData():\n")
	print(system.time({	
	minTabObj <- f.tableFromMicroData(datObj, levelObj)
	}))		
	
	cat("Time for f.tableFromMicroData():\n")
	print(system.time({	
	fullTabObj <- f.fillTable(minTabObj, levelObj, strObj)
	}))

	# set default bounds for hitas-procedure
	fullTabObj$w <- fullTabObj$Freq
	
	# cell status with possible values
	# s: potential secondary suppression
	# u: primary suppression
	# z: forced for publication	
	fullTabObj$status  <- rep("s", length(fullTabObj$strID)) 
	fullTabObj$lb  <- rep(0, length(fullTabObj$strID)) # non negative
	fullTabObj$ub <- sapply(fullTabObj$Freq, function(x) { max(2*x, 5)})
	fullTabObj$LPL <- rep(0, length(fullTabObj$strID)) # not exactly recalcable
	fullTabObj$UPL <- rep(0, length(fullTabObj$strID)) # not exactly recalcable	
	fullTabObj$SPL <- rep(1, length(fullTabObj$strID)) # not exactly recalcable 			
	
	fullTabObj$UB <- fullTabObj$ub - fullTabObj$Freq
	fullTabObj$LB <- fullTabObj$Freq - fullTabObj$lb
	
	outObj$datObj <- datObj	
	outObj$strObj <- strObj
	outObj$levelObj <- levelObj		
	outObj$fullTabObj <- fullTabObj	
	outObj
}