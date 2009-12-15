genCodesFromFile <- function(file, hierIndex="@", varname, minimal=FALSE) {
	calcDigits <- function(hier, nrLevels) {		
		# Total (always only 1 digit needed)
		nrDigits <- 1
		# level 2
		nrDigits <- c(nrDigits, nchar(as.character(length(which(hier$level==2)))))
		
		if(nrLevels > 2) {
			for(i in 2:(nrLevels-1)) {
				ss <- subset(hier, hier$level %in% c(i, i+1))
				ss$fac <- NA
				ind <- which(ss$level==i)
				for (j in 1:(length(ind)-1)) {
					ss$fac[ind[j]:(ind[j+1]-1)] <- j
				}
				spl <- split(ss, ss$fac)
				nrDigits <- c(nrDigits,	nchar(as.character((max(unlist(lapply(spl, function(x) { nrow(x)})))-1))))
			}				
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
		
		for(i in 1:nrow(hier)) {
			actLevel <- hier$level[i]
			charsActLevel <- nrDigits[actLevel] 	
			
			if(i == 1) {
				substr(hier$code[i], substrInd[[actLevel]][1], substrInd[[actLevel]][2]) <- fixedCharVec(1, charsActLevel)		
			}	
			else {
				# same level -> copy and increment
				if(hier$level[i] == hier$level[i-1]) {
					oldInd <- i-1
					hier$code[i] <- hier$code[i-1]
					oldVal <- as.integer(substr(hier$code[oldInd],substrInd[[actLevel]][1],substrInd[[actLevel]][2]))
					substr(hier$code[i], substrInd[[actLevel]][1], substrInd[[actLevel]][2]) <- fixedCharVec(oldVal+1, nrDigits[actLevel])		
				}
				
				# go back and add new hierarchy-level
				else if(hier$level[i] > hier$level[i-1]) {
					oldInd <- i-1
					hier$code[i] <- hier$code[i-1]
					oldVal <- as.integer(substr(hier$code[oldInd],substrInd[[actLevel]][1],substrInd[[actLevel]][2]))
					substr(hier$code[i], substrInd[[actLevel]][1], substrInd[[actLevel]][2]) <- fixedCharVec(oldVal+1, nrDigits[actLevel])		
				}
				
				# go back as far as necessary
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
		out
	}		
	
	hier <- read.table(file, sep=";", dec=".")	
	colnames(hier) <- c("hierarchy", varname)
	
	# calc hierarchy-level (1=TOT)
	hier$level <- as.integer(nchar(as.character(hier$hierarchy))) + 1
	
	# Nr. of Levels
	nrLevels <- length(unique(hier$hierarchy))+1 # +1 wg TOT
	
	# calculate digits
	nrDigits <- calcDigits(hier, nrLevels)
	
	# calculate standardized codes
	out <- genStandardCodes(hier, nrDigits, nrLevels)
	
	if(minimal==TRUE) 
		out <- genMinimalCodes(out)
	
	return(list(codesOrig=out[,which(colnames(out)==varname)],
				codesStandard=out$code,
				levelStructure=nrDigits))
}