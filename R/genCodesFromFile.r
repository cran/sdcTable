genCodesFromFile <- function(file, hierIndex="@", varname, minimal=FALSE) {
	calcDigits <- function(hier, hierIndex="@", nrLevels) {
		# calc split-index
		startIndex <- c(which(hier$hierarchy==hierIndex), nrow(hier))
		splitIndex <- NULL
		for (z in 1:(length(startIndex)-1)) {
			if(z != (length(startIndex)-1)) 
				splitIndex <- c(splitIndex, rep(z, length=length(startIndex[z]:(startIndex[z+1]-1))))
			else
				splitIndex <- c(splitIndex, rep(z, length=length(startIndex[z]:(startIndex[z+1]))))
		}	
		spl <- split(hier, splitIndex)
		nrDigits <- 1 # one digit needed for TOT
		nrDigits <- c(nrDigits, nchar(as.character(length(spl))))
		if(nrLevels > 2) {
			for (i in 2:(nrLevels-1)) 
				nrDigits <- c(nrDigits, nchar(as.character(max(unlist(lapply(spl, function(x) { length(which(x$level==i)) } ) )))))
		}
		nrDigits
	}	
	
	# generate standard-codes
	genStandardCodes <- function(hier, nrDigits, nrLevels) {
		# standard-digits for substrings
		genDigits <- function(nrDigits, nrLevels) {
			cs <- cumsum(nrDigits)
			substrInd <- list()
			substrInd[[1]] <- c(1,1)
			for(j in 2:nrLevels) {
				whichDigits <- (cs[(j-1)]+1):cs[j]					
				
				if(length(whichDigits)==1)
					whichDigits <- c(whichDigits,whichDigits)
				substrInd[[j]] <- whichDigits
			}
			substrInd
		}	
		
		fixedCharVec <- function(v, digits) {
			newVec <- rep(paste(rep("0",digits), collapse=""), length(v))
			for (i in 1:length(v)) {
				if(nchar(as.character(v[i])) < digits) {
					diff <- digits - nchar(as.character(v[i]))
					newVec[i] <- paste(rep("0",diff), as.character(v[i]), sep="")
				}
				else
					newVec[i] <- as.character(v[i])
			}
			newVec
		}		
		
		substrInd <- genDigits(nrDigits, nrLevels)
		hier$code <- paste(rep("0", sum(nrDigits)), collapse="")
		hier$level <- hier$level+1
		
		for(i in 1:nrow(hier)) {
			actLevel <- hier$level[i]
			charsActLevel <- nrDigits[actLevel] 	
			
			if(i == 1) {
				substr(hier$code[i], substrInd[[actLevel]][1], substrInd[[actLevel]][2]) <- fixedCharVec(1, charsActLevel)		
			}		
			else {
				# gleiches Level -> kopieren und erhöhen
				if(hier$level[i] == hier$level[i-1]) {
					oldInd <- i-1
					hier$code[i] <- hier$code[i-1]
					oldVal <- as.integer(substr(hier$code[oldInd],substrInd[[actLevel]][1],substrInd[[actLevel]][2]))
					substr(hier$code[i], substrInd[[actLevel]][1], substrInd[[actLevel]][2]) <- fixedCharVec(oldVal+1, nrDigits[actLevel])		
				}
				
				# eins zurückgehen, neue Hierarchiestufe hinzufügen
				else if(hier$level[i] > hier$level[i-1]) {
					oldInd <- i-1
					hier$code[i] <- hier$code[i-1]
					oldVal <- as.integer(substr(hier$code[oldInd],substrInd[[actLevel]][1],substrInd[[actLevel]][2]))
					substr(hier$code[i], substrInd[[actLevel]][1], substrInd[[actLevel]][2]) <- fixedCharVec(oldVal+1, nrDigits[actLevel])		
				}
				
				# soweit wie notwendig zurückgehen
				else if(hier$level[i] < hier$level[i-1]) {
					candidate <- which(hier$level==actLevel)
					oldInd <- candidate[max(which(candidate < i))]
					hier$code[i] <- hier$code[oldInd]
					oldVal <- as.integer(substr(hier$code[oldInd],substrInd[[actLevel]][1],substrInd[[actLevel]][2]))
					substr(hier$code[i], substrInd[[actLevel]][1], substrInd[[actLevel]][2]) <- fixedCharVec(oldVal+1, nrDigits[actLevel])		
				}
			}
			
		}
		hier	
	}	
	
	# remove (sub)totals
	genMinimalCodes <- function(out) {
		out$level <- nchar(as.character(out$hierarchy))
		removeInd <- NULL
		for(i in 1:(nrow(out)-1)) {
			if(out$level[i+1] > out$level[i])
				removeInd <- c(removeInd, i)
		}
		if(length(removeInd) > 0)
			out <- out[-removeInd,]
		#out	<- out[,-which(colnames(out) %in% "level")]
		#rownames(out) <- NULL
		out
	}		
	
	hier <- read.table(file, sep=";", dec=".")	
	colnames(hier) <- c("hierarchy", varname)
	
	# calc hierarchy-level (0=TOT)
	hier$level <- as.integer(nchar(as.character(hier$hierarchy)))
	
	# Nr. of Levels
	nrLevels <- length(unique(hier$hierarchy))+1 # +1 wg TOT
	
	# calculate digits
	nrDigits <- calcDigits(hier, hierIndex, nrLevels)
	
	# calculate standardized codes
	out <- genStandardCodes(hier, nrDigits, nrLevels)
	
	if(minimal==TRUE) 
		out <- genMinimalCodes(out)
	
	#rownames(out) <- NULL
	#out	<- out[,-which(colnames(out) %in% "level")]
	return(list(codesOrig=out[,which(colnames(out)==varname)],
				codesStandard=out$code,
				levelStructure=nrDigits))
}