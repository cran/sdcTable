createFullData <- function(minimalData, indexvars, l, suppVals=FALSE, suppLimit=NULL, suppZeros=NULL) {
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
	
	# fill levels that are calculated from lower hierarchies
	fillData <- function(compDat, fullDims) {
		runnInd <- FALSE
		nrIndexvars <- length(fullDims)
		vals <- compDat$val
		
		while(runnInd==FALSE) {
			for(i in 1:nrIndexvars) {			
				cat("processing dimension",i,"|",nrIndexvars," ..."); flush.console()
				tt <- Sys.time()
				for(j in (length(fullDims[[i]])):1) {				
					ss <- subset(compDat[,c(indexvars, which(colnames(compDat) == "val"))], compDat[,indexvars[i]] %in% fullDims[[i]][[j]])		
					if(length(indexvars) > 2) {	
						fac <- as.factor(apply(ss[,indexvars[-i]], 1, function(x) { paste(x, collapse="") } ))
						spl <- split(1:length(fac), fac)
					}
					else if(length(indexvars)==2){
						fac <- as.factor(as.character(ss[,indexvars[i]]))
						spl <- split(1:length(fac), fac)	
					}
					# only a single dimensional variable
					else {
						spl <- list()
						spl[[1]] <- 1:nrow(ss)
					}
					for (z in 1:length(spl)) {
						tmp <- ss[spl[[z]], ]
						ind <- which(is.na(tmp$val))
						if(length(ind) == 1)
							vals[as.numeric(rownames(tmp)[ind])] <- sum(as.integer(tmp$val), na.rm=T)	
					}				
					compDat$val <- vals	
					gc()
				}
				cat("[DONE in:", Sys.time()-tt,"]\n"); flush.console()
			}
			if(length(which(is.na(vals)))==0)
				runnInd <- TRUE
		}
		compDat
	}
	
	##########################
	if(length(indexvars) > 1)
		existingDims <- apply(minimalData[,indexvars], 2, function(x) { unique(x) } )
	else {
		existingDims <- list()
		existingDims[[1]] <- as.character(unique(minimalData[,indexvars]))
	}	
	
	newDims <- list()
	for (i in 1:length(indexvars)) {
		newDims[[i]] <- list()
		for (j in (length(l[[i]])-1):1) {
			spl <- split(minimalData[,indexvars[i]], substr(minimalData[,indexvars[i]], 1, sum(l[[i]][1:j])))
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
	
	allDims <- levels <- list()
	for (i in 1:length(indexvars)) {
		allDims[[i]] <- unique(unlist(c(existingDims[[i]], newDims[[i]])))
		levels[[i]] <- as.integer(sapply(allDims[[i]], genLevel, l[[i]]))
	}	
	
	completeData <- expand.grid(allDims)
	completeData$val <- NA
	
	if(length(indexvars) > 1) {
		completeData$fac <- apply(completeData[,indexvars], 1, function(x) { paste(x, collapse="")} )
		minimalData$fac <- apply(minimalData[,indexvars], 1, function(x) { paste(x, collapse="")} )
	}
	else {
		completeData$fac <- as.character(completeData[,indexvars])
		minimalData$fac <- as.character(minimalData[,indexvars])
	}
	
	colnames(minimalData) <- colnames(completeData)
	completeData <- rbind(completeData[which(!completeData$fac %in% minimalData$fac),], minimalData)	
	
	colNr <- min((1:ncol(minimalData))[-indexvars])
	
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
	if(length(indexvars) > 1) {
		completeData$Levs <- apply(completeData[,(length(indexvars)+2):(2*length(indexvars)+1)], 1, function(x) { paste(x, collapse="")})
	}
	else {
		completeData$Levs <- apply(completeData[,indFac:(indFac+1)], 1, function(x) { paste(x, collapse="")})
	}
	# calculate complete (splitted) level structure
	dimensions <- list()
	for (i in 1:length(indexvars))
		dimensions[[i]] <- createLevelStructure(allDims[[i]], levels[[i]], l[[i]])
	
	# fill dataset (all possible level-combinations)
	completeData <- fillData(completeData, dimensions)
	
	# keep only necessary columns
	completeData <- completeData[,which(colnames(completeData) %in% c(paste("Var", 1:length(indexvars), sep=""), "val"))]
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
	
	# do primary suppression if requested
	completeData$val <- as.integer(completeData$val)
	if(suppVals==TRUE) {
		if(is.null(suppLimit) | is.null(suppZeros)) 
			stop("both parameters suppLimit and suppZeros are necessary!")
		else {		
			if(suppZeros==TRUE) 
				completeData[which(completeData$val <= suppLimit), "geh"] <- "P"
			else 
				completeData[which(completeData$val <= suppLimit & completeData$val > 0), "geh"] <- "P"	
		}
	}	
	
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