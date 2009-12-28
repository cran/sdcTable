createFullData <- function(minimalDat=NULL, microDat=NULL, tableVars=NULL, numVar=NULL, suppRule_Freq=NULL, suppRule_P=NULL, suppRule_NK=NULL, l=l) {
	# p%-rule
	# cell unsafe if cell-total - 2 largest contributors is 
	# smaller than p% of the largest contributor
	pPercRule <- function(TotVal, cont1, cont2, p) {
		unsafe <- F
		if((TotVal - cont1 - cont2) < (p/100*cont1))
			unsafe <- T
		unsafe
	}	
	# n/k dominance rule:
	# cell unsafe if sum of the n-largest contributors is > than k% of cell total
	nkRule <- function(TotVal, sumNcont, k) {
		unsafe <- F
		if((sumNcont) > (k/100*TotVal))
			unsafe <- T
		unsafe
	}	
	# frequency rule with threshold n and indicator, if
	# zero-vals should be suppressed or not
	freqRule <- function(cellVals, n , zeroInd) {
		unsafe <- rep(F, length(cellVals))
		if(zeroInd==TRUE)
			unsafe[which(cellVals <= n)] <- T
		else
			unsafe[which(cellVals <= n & cellVals > 0)] <- T
		unsafe
	}
	
	# calculate a table from input microdata
	genTableFromMicroData <- function(microDat, tableVars, numVar, suppRule_P, suppRule_NK, suppRule_Freq) {
		# primary suppression rules
		if(!is.null(suppRule_P))
			p <- suppRule_P
		if(!is.null(suppRule_NK)) {
			nk_n <- suppRule_NK[1]
			nk_k <- suppRule_NK[2]	
		}
		if(!is.null(suppRule_Freq)) {
			freq_n <- suppRule_Freq[1]
			freq_zero <- ifelse(suppRule_Freq[2]==0, FALSE, TRUE)	
		}		
		
		indTableVars <- which(colnames(microDat) %in% tableVars)
		
		if(length(tableVars) > 1)
			microDat$fac <- apply(microDat[,indTableVars], 1, function(x) { paste(x, collapse="")})
		else
			microDat$fac <- microDat[,indTableVars]
		
		ss <- microDat[, c(indTableVars, ncol(microDat))]
		
		tt <- as.data.frame(table(microDat[,indTableVars]))
		colnames(tt)[1:(ncol(tt)-1)] <- tableVars
		
		if(length(indTableVars) > 1)
			tt$fac <- apply(tt[,1:(ncol(tt)-1)], 1, function(x) { paste(x, collapse="")})
		else
			tt$fac <- tt[,1]
		
		if(!is.null(numVar)) {	
			indNumVar <- which(colnames(microDat) %in% numVar)
			ergVal <- aggregate(microDat$numVal, by=list(microDat[,"fac"]), sum)
			colnames(ergVal) <- c("fac", "numVal")
			tt <- merge(tt, ergVal)
			
			if(!is.null(suppRule_P) | !is.null(suppRule_NK)) {
				splitMicro <- split(microDat, microDat$fac)
				totVals <- as.numeric(unlist(lapply(splitMicro, function(x) { sum(x$numVal)})))
				geh <- rep("", nrow(tt))
				
				for(i in 1:length(splitMicro)) {
					nRows <- nrow(splitMicro[[i]])
					cont <- sort(splitMicro[[i]][,indNumVar], decreasing=TRUE)
					if(!is.null(suppRule_P)){
						if(nRows >=2) {
							gehP <- pPercRule(totVals[i], cont[1], cont[2], p)
						}
					} 
					else
						gehP <- FALSE
					
					if(!is.null(suppRule_NK)) {
						if(nRows >= nk_n) 
							gehNK <- nkRule(totVals[i], sum(cont[1:nk_n]), nk_k)
					} 
					else
						gehNK <- FALSE						
					
					if(gehP==FALSE & gehNK==FALSE)
						geh[i] <- ""
					else
						geh[i] <- "P"
				}			
				tt$geh <- geh						
			}			
			
			if(!is.null(suppRule_Freq)) {
				geh <- freqRule(tt$Freq, freq_n , freq_zero)
				tt$geh <- ""
				tt$geh[which(geh==T)] <- "P"
			}			
		}
		
		tt <- tt[,-which(colnames(tt)=="fac")]	
		if(length(which(colnames(tt) == "geh"))==0)
			tt$geh <- ""
		tt	
	}	
	
	# create a complete hierarchy
	createLevelStructure <- function(allDims, levels, l) {
		output <- list()
		dat <- data.frame(dims=allDims, lev=levels)
		z <- 1
		for(i in 1:(max(dat$lev)-1)) {
			tmp <- dat[which(dat$lev==i),]
			# Top-Level
			if(nrow(tmp)==1) {
				output[[z]] <- as.character(dat[which(dat$lev %in% c(i, i+1)),"dims"])
				z <- z + 1
			}
			# we need to split
			else {
				for(j in 1:nrow(tmp)) {
					aktDim <- as.character(tmp[j,"dims"])
					aktLev <- tmp[j,"lev"]
					erg <- c(aktDim, as.character(dat[which(dat$lev==(aktLev+1) & substr(dat$dims,1,sum(l[1:aktLev])) == substr(aktDim,1,sum(l[1:aktLev]))), "dims"]))	
					if(length(erg)>1) {
						output[[z]] <- erg	
						z <- z + 1
					}
				}
			}
		}
		output
	}
	
	# calculate Levels
	genLevel <- function(level, dimStructure) {
		cums <- cumsum(dimStructure)
		lenStructure <- length(dimStructure)
		if(as.integer(substr(level, 1,cums[lenStructure]))==0)
			out <- 1
		else if(as.integer(substr(level, cums[lenStructure-1]+1, cums[lenStructure]))!=0)
			out <- lenStructure
		else {
			for(i in (2:(lenStructure-1))) {
				if(as.integer(substr(level, cums[i-1]+1, cums[i]))!=0 & as.integer(substr(level, cums[i]+1, cums[lenStructure]))==0)						
					out <- i
			}				
		}
		out
	}	
	
	# check, ob alle Merkmalskombinationen in der Tablle vorkommen
	checkMinimalData <- function(minimalData, existingDims, indexvars) {
		gridExistingDims <- expand.grid(existingDims)
		colnames(gridExistingDims) <- colnames(minimalData)[indexvars]
		if(nrow(minimalData)!=nrow(gridExistingDims)) {
			if(length(indexvars) > 1) 
				gridExistingDims$fac <- apply(gridExistingDims[,indexvars], 1, function(x) { paste(x, collapse="")} )
			else 
				gridExistingDims$fac <- as.character(gridExistingDims[,indexvars])
			
			gridExistingDims <- merge(gridExistingDims, minimalData, by="fac", all.x=T)
			
			indRem <- which(substr(colnames(gridExistingDims),nchar(colnames(gridExistingDims)),nchar(colnames(gridExistingDims)))=="y")
			gridExistingDims <- gridExistingDims[,-indRem]
			indFac <- which(colnames(gridExistingDims)=="fac")
			indNonFac <- setdiff(1:ncol(gridExistingDims), indFac)
			gridExistingDims <- gridExistingDims[,c(indNonFac, indFac)]
			colnames(gridExistingDims) <- colnames(minimalData)
			
			indNA <- which(is.na(gridExistingDims$Freq))	
			gridExistingDims$Freq[indNA] <- 0
			gridExistingDims$numVal[indNA] <- 0			
			gridExistingDims$geh[indNA] <- ""	
			minimalData <- gridExistingDims
		}		
		minimalData		
	}		
	
	# fill levels that are calculated from lower hierarchies + 
	# if necessary, primary suppression
	fillData <- function(compDat, fullDims, exDims, suppRule_Freq, suppRule_P, suppRule_NK, microDat, indexvars) {
		calcMicroDataforTotal <- function(existingDims, actDims, microData, indexvars) {
			indexVals <- list()
			str <- paste("microDat[which(", sep="")
			for (i in 1:length(actDims)) {				
				x <- unlist(strsplit(actDims[i],""))
				if(sum(as.integer(x))==0)
					indexVals[[i]] <- existingDims[[i]]
				else {
					ind <- max(which(x=="0"))-1
					if(ind == 0) {
						indexVals[[i]] <- actDims[i]
					}
					else if(ind > 1)
						indexVals[[i]] <- existingDims[[i]][which(substr(existingDims[[i]],1,ind) %in% substr(actDims[i],1,ind))]
					else 
						indexVals[[i]] <- existingDims[[i]]					
				}
				
				if(i==1)
					str <- paste(str, "microDat[,indexvars[",i,"]] %in% indexVals[[",i,"]]",sep="")
				else
					str <- paste(str, "& microDat[,indexvars[",i,"]] %in% indexVals[[",i,"]]",sep="")
				
			}
			str <- paste(str,"),]")
			micro <- eval(parse(text=str))
			micro
		}
		
		# Parameter for primary suppression rules
		if(!is.null(suppRule_P))
			p <- suppRule_P
		if(!is.null(suppRule_NK)) {
			nk_n <- suppRule_NK[1]
			nk_k <- suppRule_NK[2]
		}	
		if(!is.null(suppRule_Freq)) {
			freq_n <- suppRule_Freq[1]
			freq_zero <- ifelse(suppRule_Freq[2]==0, FALSE, TRUE)	
		}				
		
		numVal <- FALSE
		if(length(which(colnames(compDat)=="numVal"))> 0)
			numVal <- TRUE
		
		compDat <- as.matrix(compDat)
		runnInd <- FALSE
		nrIndexvars <- length(fullDims)
		vals <- as.numeric(compDat[,"Freq"])
		if(length(which(colnames(compDat)=="numVal"))> 0)
			valsNum <- as.numeric(compDat[,"numVal"])
		gehVec <- compDat[,"geh"]
		while(runnInd==FALSE) {
			for(i in 1:nrIndexvars) {			
				cat("processing dimension",i,"|",nrIndexvars," ..."); flush.console()
				tt <- Sys.time()
				for(j in (length(fullDims[[i]])):1) {	
					indRows <- which(compDat[,indexvars[i]] %in% fullDims[[i]][[j]])
					ss <- compDat[indRows,c(indexvars, which(colnames(compDat) %in% c("Freq","numVal", "fac")))]
					if(length(indexvars) > 2) {	
						fac <- as.factor(apply(ss[,indexvars[-i]], 1, function(x) { paste(x, collapse="") } ))
						spl <- split(1:length(fac), fac)
					}
					else if(length(indexvars)==2){
						fac <- as.factor(as.character(ss[,indexvars[-i]]))
						spl <- split(1:length(fac), fac)	
					}
					# only a single dimensional variable
					else {
						spl <- list()
						spl[[1]] <- 1:nrow(ss)
					}
					
					for (z in 1:length(spl)) {
						tmp <- ss[spl[[z]], ]						
						ind <- which(is.na(as.numeric(tmp[,"Freq"])))
						if(length(ind) == 1) {
							totFreq <- sum(as.integer(tmp[,"Freq"]), na.rm=T)
							if(numVal==TRUE)
								totVal <- sum(as.integer(tmp[,"numVal"]), na.rm=T)
							rowInd <- indRows[spl[[z]][ind]]
							vals[rowInd] <- totFreq
							if(numVal==TRUE)
								valsNum[rowInd] <- totVal
							nRows <- nrow(tmp) - 1
							
							# only useful if cell is populated
							if(sum(as.integer(tmp[!is.na(tmp[,"Freq"]),"Freq"]))>1 & numVal==TRUE) {
								if(!is.null(suppRule_P) | !is.null(suppRule_NK)) {
									subMicroDat <- calcMicroDataforTotal(exDims, tmp[which(is.na(tmp[,"Freq"])),indexvars], microDat, indexvars)
									contributions <- sort(subMicroDat[,"numVal"], decreasing=T)
									if(!is.null(suppRule_P)) {
										if(nRows >= 2) {
											gehP <- pPercRule(totVal, contributions[1], contributions[2], p)						
										}	
									}
									else 
										gehP <- FALSE	
									
									if(!is.null(suppRule_NK)) {
										if(nRows >= nk_n) 
											gehNK <- nkRule(totVal, sum(contributions[1:nk_n]), nk_k)					
									}
									else 
										gehNK <- FALSE	
									
									if(gehP==FALSE & gehNK==FALSE)
										gehVec[rowInd] <- ""
									else
										gehVec[rowInd] <- "P"
								}	
							}					
						}
					}	
					compDat[,"Freq"] <- as.character(vals)
					if(numVal==TRUE)
						compDat[,"numVal"] <- as.character(valsNum)
					compDat[,"geh"] <- as.character(gehVec)
				}
				cat("[DONE in:", round(Sys.time()-tt,2),"]\n"); flush.console()
			}
			if(length(which(is.na(vals)))==0) 
				runnInd <- TRUE
		}
		
		if(!is.null(suppRule_Freq)) {
			geh <- freqRule(as.integer(compDat[,"Freq"]), freq_n , freq_zero)
			compDat[,"geh"][which(geh==T)] <- "P"
		}		
		
		compDat <- as.data.frame(compDat)
		compDat
	}	
	
	### start function ###	
	if(is.null(minimalDat) & is.null(microDat))
		stop("You need to specify either a table (minimalDat) or microData as an input object!\n")
	else if(!is.null(minimalDat) & !is.null(microDat)) {
		stop("You cannot spefify both a table and microData as input objects.\n")
	}
	else {
		# we can't calculate dominance rules
		if(!is.null(minimalDat)) {			
			suppRule_P <- NULL
			suppRule_NK <- NULL		
		}
	}	
	
	if(is.null(suppRule_P) & is.null(suppRule_NK) & is.null(suppRule_Freq)) 
		cat("Warning: no rules for primary suppression have been specified!\n")
	
	# only calculate table from microdata if microdata have
	# been speficied as input
	if(!is.null(microDat)) 
		minimalDat <- genTableFromMicroData(microDat, tableVars=tableVars, numVar=numVar, suppRule_P=suppRule_P, suppRule_NK=suppRule_NK, suppRule_Freq=suppRule_Freq) 
	else {
		if(length(which(colnames(minimalDat) =="geh"))==0)
			minimalDat$geh <- ""
	}
	
	if(!is.null(numVar))
		colnames(minimalDat)[colnames(minimalDat)==numVar] <- "numVal"
	
	indexvars <- which(colnames(minimalDat) %in% tableVars)
	
	# get existing characteristics of dimensional variables
	if(length(indexvars) > 1) {
		ex <- apply(minimalDat[,indexvars], 2, function(x) { unique(x) } )
		if(is.list(ex))
			existingDims <- ex
		else {
			existingDims <- list()
			for(i in 1:length(indexvars))
				existingDims[[i]] <- ex[,i]			
		}
	}
	else {
		existingDims <- list()
		existingDims[[1]] <- as.character(unique(minimalDat[,indexvars]))
	}		
	
	# calculate all possible characteristics (=sub|totals) 
	# for each dimensional variable
	newDims <- list()
	for (i in 1:length(indexvars)) {
		newDims[[i]] <- list()
		for (j in (length(l[[i]])-1):1) {
			spl <- split(minimalDat[,indexvars[i]], substr(minimalDat[,indexvars[i]], 1, sum(l[[i]][1:j])))
			spl <- lapply(spl, function(x) { as.character(unique(x))[1] } )
			for(z in 1:length(spl)) {
				# we calculate the upper limit and update the levels
				upperHier <- spl[[z]]
				from <- sum(l[[i]][1:j]) +1
				to <- nchar(upperHier)
				substr(upperHier, from, to) <- paste(rep("0", (to - from + 1)), collapse="")
				if(!upperHier %in% as.character(existingDims[[i]])) 					
					newDims[[i]] <- c(newDims[[i]], upperHier)			
			}
		}	
		# add new hierarchies
		if(length(newDims[[i]]) > 0)
			newDims[[i]] <- unlist(newDims[[i]])			
	}
	
	# combine existing and possible new characteristics
	allDims <- levels <- list()
	for (i in 1:length(indexvars)) {
		allDims[[i]] <- unique(unlist(c(existingDims[[i]], newDims[[i]])))
		levels[[i]] <- as.integer(sapply(allDims[[i]], genLevel, l[[i]]))
	}	
	
	# generate complete Dataset
	completeData <- expand.grid(allDims)
	completeData$Freq <- NA
	
	if(!is.null(numVar))
		completeData$numVal <- NA
	
	completeData$geh <- ""	
	if(length(indexvars) > 1) {
		completeData$fac <- apply(completeData[,indexvars], 1, function(x) { paste(x, collapse="")} )
		minimalDat$fac <- apply(minimalDat[,indexvars], 1, function(x) { paste(x, collapse="")} )
	}
	else {
		completeData$fac <- as.character(completeData[,indexvars])
		minimalDat$fac <- as.character(minimalDat[,indexvars])
	}	
	
	# check if minimalDat contains all possible combinations of minimal spanning variables
	minimalDat <- checkMinimalData(minimalDat, existingDims, indexvars)
	
	colnames(minimalDat) <- colnames(completeData)
	completeData <- rbind(completeData[which(!completeData$fac %in% minimalDat$fac),], minimalDat)	
	
	# add levels
	tmpDat <- list()	
	for (i in 1:length(indexvars)) {
		tmpDat[[i]] <- data.frame(allDims[[i]], levels[[i]])
		colnames(tmpDat[[i]]) <- c(paste("Var",i,sep=""), paste("Level",i,sep=""))
		completeData <- merge(completeData, tmpDat[[i]])
	}
	completeData <- completeData[,c(rev(which(substr(colnames(completeData),1,3) == "Var")), (length(indexvars)+1):ncol(completeData))]	
	rm(tmpDat)	
	
	indFac <- which(colnames(completeData)=="fac")
	if(length(indexvars) > 1) 
		completeData$Levs <- apply(completeData[,(length(indexvars)+2):(2*length(indexvars)+1)], 1, function(x) { paste(x, collapse="")})
	else 
		completeData$Levs <- apply(completeData[,indFac:(indFac+1)], 1, function(x) { paste(x, collapse="")})
	
	# calculate complete (splitted) level structure
	dimensions <- list()
	for (i in 1:length(indexvars))
		dimensions[[i]] <- createLevelStructure(allDims[[i]], levels[[i]], l[[i]])
	
	# fill dataset (all possible level-combinations)
	completeData <- fillData(compDat=completeData, fullDims=dimensions, exDims=existingDims, suppRule_Freq=suppRule_Freq, suppRule_P=suppRule_P, suppRule_NK=suppRule_NK, microDat=microDat, indexvars=indexvars)
	
	# keep only necessary columns
	if(is.null(numVar)) {
		completeData <- completeData[,which(colnames(completeData) %in% c(paste("Var", 1:length(indexvars), sep=""), "Freq","geh"))]
		colnames(completeData)[colnames(completeData)=="Freq"] <- "val"
	}
	else {
		completeData <- completeData[,which(colnames(completeData) %in% c(paste("Var", 1:length(indexvars), sep=""), "numVal","geh"))]
		colnames(completeData)[colnames(completeData)=="numVal"] <- "val"
		
	}
	if(length(which(colnames(completeData)=="geh"))==0)
		completeData$geh <- ""
	
	# sort by index-variables
	str <- paste("completeData <- completeData[order(")
	for (i in 1:length(indexvars)) {
		if(i < length(indexvars))
			str <- paste(str, "as.character(completeData$", colnames(completeData)[indexvars[i]], "), ", sep="")
		else 
			str <- paste(str, "as.character(completeData$", colnames(completeData)[indexvars[i]], ")),]", sep="")
	}
	eval(parse(text=str))
	rownames(completeData) <- 1:nrow(completeData)
	completeData$geh <- as.character(completeData$geh)
	completeData$val <- as.numeric(as.character(completeData$val))
	
	# create object of class "fullData
	fullData <- list()
	fullData$data <- completeData
	fullData$dims <- l
	fullData$indexvars <- indexvars
	fullData$supps2check <- which(completeData$geh != "")
	fullData$numberindexvars <- length(indexvars)
	fullData$allDims <- allDims
	fullData$dimensions <- dimensions
	class(fullData) <- "fullData"
	return(fullData)
}