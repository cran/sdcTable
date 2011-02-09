###################################################################
### calculate all necessary and possible meta-information		###
### given any valid input-object (file or dataframe) in advance ###
###################################################################
calcDimInfos <- function(inputDat, file=NULL, dataframe=NULL, vName) {
	# add the variable-name of the hierarchy and its position in the 
	# data-set to to levelInfo -> neccessary!
	# TODO: is this used multiple times or only in hierInfo()?
	f.setPosition <- function(levelInfo, inputDat, vName) {
		if(class(levelInfo) != "levelInfo")
			stop("please use an appropriate input object derived from f.hierInfo()\n")
		
		posIndex <- which( colnames(inputDat) == vName )
		if( length(posIndex) != 1 ) 
			stop("please check if 'vname' exists in 'inputDat'!\n")
		else {
			levelInfo$varName <- vName		
			levelInfo$posIndex <- posIndex		
		}		
		levelInfo
	}	
	
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
	
	f.cleanStructure <- function(dataframe, hierIndex, out) {
		hier <- dataframe
		colnames(hier) <- c("hierarchy", vName)	
		hier[,1] <- as.character(hier[,1])
		hier[,2] <- as.character(hier[,2])
		if ( length(out$levelsRemoveOrig) > 0 ) {
			for ( i in 1:length(out$levelsRemoveOrig) )	{
				# remove level from hier
				index <- which(hier[,2] == out$levelsRemoveOrig[i])
				if ( nchar(hier[index,1]) > nchar(hier[index+1,1]) ) 
					hier <- hier[-index,]
				else {
					x1 <- which(nchar(hier[,1]) == nchar(hier[index+1,1]))
					x1 <- x1[x1>index]
					tmpResult <- diff(x1)
					levChangeInd <- (index+1):(index+1+min(which(tmpResult!=1))-1)
					hier[levChangeInd,1] <- substr(hier[levChangeInd,1],2,nchar(hier[index+1,1]))
					hier <- hier[-index,]
				}				
			}
		}
		
		#correct levels:
		if ( any(diff(nchar(hier[,1]))>1) ) {
			startInd <- which(diff(nchar(hier[,1]))>1)+1
			for ( i in 1:length(startInd) ) {
				curLev <- nchar(hier[startInd[i],1])
				runInd <- TRUE
				x <- startInd[i]
				while (runInd==TRUE) {
					if ( nchar(hier[x,1]) == curLev )
						hier[x,1]	<- substr(hier[x,1], 2, nchar(hier[x,1]))
					else
						runInd <- FALSE
					x <- x+1					
				}				
			}
			
		}		
		hier	
	}
	
	f.simplify <- function(out) {
		changes <- FALSE	
		dims <- out$dimensions
		levelsRemove <- levelsRemoveUp <- levelsRemoveOrig <- levelsRemoveOrigUp <- NULL
		indexRemove <- NULL
		# levels with length 2
		lev.len2 <- which(unlist(lapply(dims, length))==2)
		if( length(lev.len2) > 0 ) {
			changes <- TRUE
			nonNecc <- NULL
			for ( i in length(lev.len2):1 ) {
				up <- dims[[lev.len2[i]]][1]
				down <- dims[[lev.len2[i]]][2]	
				
				# is the current level available in the original levels?
				indUp <- which(out$codesStandard == up)
				indDown <- which(out$codesStandard == down)
				
				## both levels are available in the original levels
				if( length(indUp) == 1 & length(indDown) == 1 ) {
					# save and remove the lower level
					levelsRemove <- c(levelsRemove, as.character(out$codesStandard)[indDown])
					levelsRemoveUp <- c(levelsRemoveUp, as.character(out$codesStandard)[indUp])
					levelsRemoveOrig <- c(levelsRemoveOrig, as.character(out$codesOrig)[indDown])
					levelsRemoveOrigUp <- c(levelsRemoveOrigUp, as.character(out$codesOrig)[indUp])
					
					indexRemove <- c(indexRemove, indDown)
					
					### aus testweise
				}
				## only one level is available in the original levels
				else {
					# if the lower level is listed: change the levels!
					if( length(indDown) == 1 )
						out$levelsOrig[indDown] <- 	out$levelsOrig[indDown]-1
				}
			}
		}
		if ( length(indexRemove) > 0 ) {
			out$codesOrig <- out$codesOrig[-indexRemove]
			out$codesStandard <- out$codesStandard[-indexRemove]
			out$levelsOrig <- out$levelsOrig[-indexRemove]				
		}	
		
		outObj <- list(
				out=out, 
				levelsRemove=levelsRemove, 
				levelsRemoveUp=levelsRemoveUp,
				levelsRemoveOrig=levelsRemoveOrig, 
				levelsRemoveOrigUp=levelsRemoveOrigUp, 
				changes=changes)
		outObj
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
	dataframe <- as.data.frame(hh)
	colnames(dataframe) <-  c("hierarchy", vName)
	file <- NULL
	
	out <- f.calcInfo(dataframe, hierIndex, vName)
	levelsRemove <- levelsRemoveUp <- levelsRemoveOrig <- levelsRemoveOrigUp <- NULL
	changes <- TRUE
	while( changes==TRUE ) {
		outS <- f.simplify(out)	
		if( outS$changes == TRUE ) {
			df <- data.frame(x=sapply(outS$out$levelsOrig, function(x) { paste(rep("@",x), collapse="") } ), y=outS$out$codesOrig)
			out <- f.calcInfo(dataframe=df, hierIndex, vName)	
			levelsRemove <- c(levelsRemove, outS$levelsRemove)
			levelsRemoveUp <- c(levelsRemoveUp, outS$levelsRemoveUp)
			levelsRemoveOrig <- c(levelsRemoveOrig, outS$levelsRemoveOrig)	
			levelsRemoveOrigUp <- c(levelsRemoveOrigUp, outS$levelsRemoveOrigUp)			
		}
		else 
			changes <- FALSE
	}
	
	if ( !is.null(outS) ) {
		out <- outS$out	
		out$levelsRemoveOrig <- levelsRemoveOrig
		out$levelsRemoveOrigUp <- levelsRemoveOrigUp	
	}
	
	df <- f.cleanStructure(dataframe, hierIndex, out)
	out <- f.calcInfo(dataframe=df, hierIndex, vName)
	
	if( length(levelsRemoveOrig) > 0 ) {
		ls <- out$levelStructure
		cs <- cumsum(ls)
		
		levRem <- levUp <- codeRem <- codeUp <- NULL
		#ordering <- order(levelsRemoveOrig)
		ordering <- order(levelsRemove)
	
		levelsRemoveOrig <- levelsRemoveOrig[ordering]
		levelsRemoveOrigUp <- levelsRemoveOrigUp[ordering]
		
		for ( i in 1:length(levelsRemoveOrig) ) {
			xx <- list()
			curLevel <- levelsRemoveOrigUp[i]			
			index <- which(out$codesOrig == curLevel)
			
			# check if upper level is removed too
			if ( length(index) == 0 ) {
				index2 <- which(levRem == curLevel)
				if (length(index2) == 0)
					stop("error in calcHierarchyInfos(), i=",i,"!\n")
				codeDown <- codeRem[index2]
				digit <- max(which(unlist(strsplit(codeDown, "")) != "0")) + 1
				from <- cs[max(which(cs < digit))] + 1
				to <- cs[min(which(cs >= from))]
				if ( to - from == 0 )
					substr(codeDown, from, to) <- "1"
				else
					substr(codeDown, from, to) <- paste(paste(rep("0", to - from), collapse=""), "1", sep="")
				xx[[1]] <- levelsRemoveOrig[i]
				xx[[2]] <- levelsRemoveOrigUp[i]
				xx[[3]] <- codeUp[index2]
				xx[[4]] <- codeDown				
			} 
			else {
				codeDown <- out$codesStandard[index]
				digit <- max(which(unlist(strsplit(codeDown, "")) != "0")) + 1
				from <- cs[max(which(cs < digit))] + 1
				to <- cs[min(which(cs >= from))]
				if ( to - from == 0)
					substr(codeDown, from, to) <- "1"
				else
					substr(codeDown, from, to) <- paste(paste(rep("0", to - from), collapse=""), "1", sep="")
				xx[[1]] <- levelsRemoveOrig[i]
				xx[[2]] <- levelsRemoveOrigUp[i]
				xx[[3]] <- out$codesStandard[index]
				xx[[4]] <- codeDown				
			}	
			levRem <- c(levRem, xx[[1]])
			levUp <- c(levUp, xx[[2]])
			codeUp <- c(codeUp, xx[[3]])
			codeRem <- c(codeRem, xx[[4]])				
		}			
		out$levelsRemoveOrig <- levRem
		out$levelsRemoveOrigUp <- levUp
		out$codeRemoveOrig <- codeRem
		out$codeRemoveOrigUp <- codeUp
	}	
	class(out) <- "levelInfo"	
	out <- f.setPosition(out, inputDat, vName)
	out
}