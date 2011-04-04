###################################################################
### calculate all necessary and possible meta-information		###
### given any valid input-object (file or dataframe) in advance ###
###################################################################
calcDimInfos <- function(inputDat, file=NULL, dataframe=NULL, vName) {
	f.calcInfo <- function(dataframe, hierIndex, vName) {
		genLevel <- function(level, dimStructure) {
			cums <- cumsum(dimStructure)
			lenStructure <- length(dimStructure)
			if( as.integer(substr(level, 1, cums[lenStructure])) == 0 )
				out <- 1
			else if( as.integer(substr(level, cums[lenStructure-1]+1, cums[lenStructure])) != 0 )
				out <- lenStructure
			else {
				for( i in (2:(lenStructure-1)) ) {
					if( as.integer(substr(level, cums[i-1]+1, cums[i]))!=0 & as.integer(substr(level, cums[i]+1, cums[lenStructure])) == 0 )						
						out <- i
				}				
			}
			out
		}	
		
		hier <- dataframe
		colnames(hier) <- c("hierarchy", vName)	
		
		# 2) calculate the levels and the the number of levels
		hier$level <- as.integer(nchar(as.character(hier$hierarchy))) 	
		nrLevels <- length(unique(hier$hierarchy)) 
		
		# 2) calculate necessary digits to represent this hierarchy
		nrDigits <- c(1, nchar(as.character(length(which(hier$level==2))))) # level 2
		if( nrLevels > 2 ) {
			for( i in 2:(nrLevels-1) ) {
				ss <- subset(hier, hier$level %in% c(i, i+1))
				ss$fac <- NA
				ind <- which(ss$level==i)
				for ( j in 1:(length(ind)) ) {
					if ( j != length(ind) ) 
						ss$fac[ind[j]:(ind[j+1]-1)] <- j		
					else 
						ss$fac[ind[j]:length(ss$fac)] <- j	
				} 								
				spl <- split(ss, ss$fac)
				nrDigits <- c(nrDigits,	nchar(as.character((max(unlist(lapply(spl, function(x) { nrow(x)})))-1))))
			}				
		}	
		
		# 3) calculate standard-codes
		# calculate position of levels in standard codec (substrInd)
		cs <- cumsum(nrDigits)
		substrInd <- list()
		substrInd[[1]] <- c(1,1)
		for( j in 2:nrLevels ) {
			whichDigits <- (cs[(j-1)]+1):cs[j]					
			if( length(whichDigits) == 1 )
				whichDigits <- c(whichDigits,whichDigits)
			substrInd[[j]] <- whichDigits
		}
		
		hier$code <- paste(rep("0", sum(nrDigits)), collapse="")
		# calc the standard codes
		for( i in 1:nrow(hier) ) {
			actLevel <- hier$level[i]
			charsActLevel <- nrDigits[actLevel] 	
			
			if( i == 1 ) {
				if ( nchar(hier$hierarchy[1]) != 1 )	
					substr(hier$code[i], substrInd[[actLevel]][1], substrInd[[actLevel]][length(substrInd[[actLevel]])]) <- sprintf(paste("%0",charsActLevel,"d",sep=""),1)
			}				
			else {
				if( hier$level[i] >= hier$level[i-1] ) {
					oldInd <- i-1
					hier$code[i] <- hier$code[i-1]
					oldVal <- as.integer(substr(hier$code[oldInd],substrInd[[actLevel]][1],substrInd[[actLevel]][length(substrInd[[actLevel]])]))
					substr(hier$code[i], substrInd[[actLevel]][1], substrInd[[actLevel]][length(substrInd[[actLevel]])]) <- sprintf(paste("%0",charsActLevel,"d",sep=""),oldVal+1)
				}	
				# go back as far as necessary
				else if( hier$level[i] < hier$level[i-1] ) {
					candidate <- which(hier$level==actLevel)
					oldInd <- candidate[max(which(candidate < i))]
					hier$code[i] <- hier$code[oldInd]
					oldVal <- as.integer(substr(hier$code[oldInd],substrInd[[actLevel]][1],substrInd[[actLevel]][length(substrInd[[actLevel]])]))
					substr(hier$code[i], substrInd[[actLevel]][1], substrInd[[actLevel]][length(substrInd[[actLevel]])]) <- sprintf(paste("%0",charsActLevel,"d",sep=""),oldVal+1)
				}
			}			
		}
		
		# calculate if a given level is neccessary
		level <- hier$level
		removeInd <- NULL
		for( i in 1:(nrow(hier)-1) ) {
			if( level[i+1] > level[i] )
				removeInd <- c(removeInd, i)
		}
		hier$codesMinimal <- TRUE
		if( length(removeInd) > 0 )
			hier$codesMinimal[removeInd] <- FALSE
		
		codesOrig <- as.character(hier[,colnames(hier)==vName])
		codesStandard <- as.character(hier$code)
		levelsOrig <- hier$level
		codesMinimal <- as.character(hier$codesMinimal)
		levelStructure <- nrDigits		
		
		# 4) calculate additional information
		minInd <- codesStandard[codesMinimal==TRUE]
		
		# calculate all possible characteristics (=sub|totals) 
		# for each dimensional variable
		newDims <- NULL
		for ( j in (length(nrDigits)-1):1 ) {
			spl <- split(minInd, substr(minInd, 1, sum(nrDigits[1:j])))
			spl <- lapply(spl, function(x) { as.character(unique(x))[1] } )
			for( z in 1:length(spl) ) {
				# we calculate the upper limit and update the levels
				upperHier <- spl[[z]]
				from <- sum(nrDigits[1:j]) + 1
				to <- nchar(upperHier)
				substr(upperHier, from, to) <- paste(rep("0", (to - from + 1)), collapse="")
				if( !upperHier %in% as.character(codesStandard) ) 					
					newDims <- append(newDims, upperHier)			
			}
		}	
		
		# combine existing and possible new characteristics
		allDims <- unique(unlist(c(codesStandard, newDims)))
		levels <- as.integer(sapply(allDims, genLevel, nrDigits))
		
		dimensions <- list()
		dat <- data.frame(dims=allDims, lev=levels)
		z <- 1
		for( i in 1:(max(dat$lev)-1) ) {
			tmp <- dat[which(dat$lev==i),]
			# Top-Level
			if( nrow(tmp) == 1 ) {
				dimensions[[z]] <- as.character(dat[which(dat$lev %in% c(i, i+1)),"dims"])
				z <- z + 1
			}
			# we need to split
			else {
				for( j in 1:nrow(tmp) ) {
					aktDim <- as.character(tmp[j,"dims"])
					aktLev <- tmp[j,"lev"]
					erg <- c(aktDim, as.character(dat[which(dat$lev==(aktLev+1) & substr(dat$dims,1,sum(nrDigits[1:aktLev])) == substr(aktDim,1,sum(nrDigits[1:aktLev]))), "dims"]))	
					if( length(erg) > 1 ) {
						dimensions[[z]] <- erg	
						z <- z + 1
					}
				}
			}
		}
		dimensions <- lapply(dimensions, sort)
		
		# recalculate the neccessary (TRUE) and non-neccessary levels
		# based on all possible levels (out$allDims)
		notUsed <- sort(unique(unlist(lapply(dimensions, function(x) x[1]))))
		codesMinimal <- rep(TRUE, length(allDims))
		codesMinimal[allDims %in% notUsed] <- FALSE		
		
		out <- list(
				codesOrig=codesOrig,
				codesStandard=codesStandard,
				levelsOrig=levels,
				levelStructure=nrDigits,
				#allDims=allDims,
				#levels=levels,
				dimensions=dimensions,
				codesMinimal=codesMinimal)		
		out		
	}	
	
	hierIndex <- "@"
	if( is.null(file) & is.null(dataframe) )
		stop("Please specify either a file or a data.frame with the hierachy-codes!\n")
	
	if( !is.null(file) ) {
		hh <- read.table(file, sep=";", dec=".", colClasses="character")
		if( ncol(hh) > 2 )
			stop("Please use a correct input file! (Hint: 2 columns only!)\n")
	}
	if( !is.null(dataframe) ) {
		if( ncol(dataframe) > 2 )
			stop("Please use a correct input dataframe! (Hint: 2 columns only!)\n")
		hh <- dataframe
	}		
	hh <- as.matrix(hh)		
	hh[,1] <- as.character(hh[,1])
	hh[,2] <- as.character(hh[,2])			
	
	if ( substr(hh[1,1],1,1) != "@" )
		stop("please check your input data. The character specifiying level-information must be '@'!\n")
	
	# if grand total (eg. 00000x) is not specified! -> "@" is the grand total
	if ( nchar(as.character(hh[1,1])) > 1 ) 
		hh <- rbind(c(hierIndex,"TOT"), hh)
	
	#hh <- subset(hh, substr(hh[,2],1,1)=="D")
 	#hh[,1] <- sapply(hh[,1], function(x) { substr(x, 2, nchar(x)) } )
	dataframe <- as.data.frame(hh)
	colnames(dataframe) <-  c("hierarchy", vName)
	dataframe[,1] <- as.character(dataframe[,1])
	file <- NULL
	
	### calc entire level-structure ###
	infoComplete <- f.calcInfo(dataframe, hierIndex, vName)
	
	### search for duplicates ####
	dimLen <- sapply(infoComplete$dimensions, length)
	dups <- dupsUp <- NULL
	removeInd <- NULL
	if ( any(dimLen == 2 ) ) {
		index <- which(dimLen==2)	
		for ( i in 1:length(index)) {
			indexInOrig1 <- match(infoComplete$dimensions[[index[i]]][2], infoComplete$codesStandard)
			levDiff <- setdiff(which(infoComplete$levels < infoComplete$levels[indexInOrig1]), 1:indexInOrig1)
			if ( length(levDiff) > 0 ) 
				indexInOrig2 <- min(levDiff)				
			else 
				indexInOrig2 <- nrow(dataframe)+1
			
			# move one level up
			if ( indexInOrig2 - indexInOrig1 > 1 ) {
				ind <- (indexInOrig1+1):(indexInOrig2-1)
				dataframe[ind,1] <- substr(dataframe[ind,1], 2, nchar(dataframe[ind,1]))
			}		
	
			# add info
			dups <- c(dups, infoComplete$codesOrig[indexInOrig1])
			dupsUp <- c(dupsUp, infoComplete$codesOrig[indexInOrig1-1])
			removeInd <- c(removeInd, indexInOrig1)
		}
		dataframe <- dataframe[-removeInd,]
	}
	
	info <- f.calcInfo(dataframe, hierIndex, vName)
	info$dups <- dups
	info$dupsUp <- dupsUp
	
	# set position
	posIndex <- which( colnames(inputDat) == vName )
	if( length(posIndex) != 1 ) 
		stop("please check if 'vname' exists in 'inputDat'!\n")
	else {
		info$varName <- vName		
		info$posIndex <- posIndex		
	}		
	info
}
