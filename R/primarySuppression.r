###########################################
### doing primary suppression according ###
### to defined rules on a meta-object   ###
###########################################
primarySuppression <- function(outObj, suppRule_Freq=NULL, suppRule_P=NULL, suppRule_NK=NULL) {
	# p%-rule
	# cell unsafe if cell-total - 2 largest contributors is 
	# smaller than p% of the largest contributor
	f.suppRule.pPerc <- function(celltot, cont1, cont2, p) {
		unsafe <- FALSE
		if( (celltot - cont1 - cont2) < (p/100*cont1) )
			unsafe <- TRUE
		unsafe
	}	
	# n/k dominance rule:
	# cell unsafe if sum of the n-largest contributors is 
	# larger than k% of cell total
	f.suppRule.nk <- function(celltot, sumNcont, k) {
		unsafe <- FALSE
		if( (sumNcont) > (k/100*celltot) )
			unsafe <- TRUE
		unsafe
	}	
	# frequency rule with threshold n and indicator, if
	# zero-vals should be suppressed or not
	f.suppRule.freq <- function(cellVals, n, zeroInd) {
		unsafe <- rep(F, length(cellVals))
		if( zeroInd==TRUE )
			unsafe[which(cellVals <= n)] <- T
		else
			unsafe[which(cellVals <= n & cellVals > 0)] <- T
		unsafe
	}		
	
	fullTabObj <- outObj$fullTabObj
	datObj <- outObj$datObj
	strObj <- outObj$strObj	
	levelObj <- outObj$levelObj
	suppsP <- suppsNK <- suppsFreq <- rep(FALSE, length(fullTabObj$Freq))		
	
	# may be used with or without microdata
	if( !is.null(suppRule_Freq) ) {
		freq_n <- suppRule_Freq[1]
		freq_zero <- ifelse(suppRule_Freq[2]==0, FALSE, TRUE)	
		suppsFreq <- f.suppRule.freq(as.integer(fullTabObj$Freq), freq_n , freq_zero)
	}		
	if( datObj$microDat == FALSE & any(!is.null(c(suppRule_P, suppRule_NK))) ) { 
		cat("--> Info: microdata are not available. We can't use any dominance rules!\n")
		suppRule_P 	<- NULL
		suppRule_NK <- NULL		
	}		
	
	switch <- TRUE
	if ( any(is.na(fullTabObj$numVal)) & any(!is.null(c(suppRule_P, suppRule_NK))) ) {
		cat("Info: Dominance rules can't be used because no numerical variable is available!\n")
		suppRule_P 	<- NULL
		suppRule_NK <- NULL
		switch <- FALSE
	}	
	
	if ( any(!is.null(c(suppRule_P, suppRule_NK))) & switch==TRUE ) {
		if( !is.null(suppRule_P) )
			p <- suppRule_P
		if( !is.null(suppRule_NK) ) {
			nk_n <- suppRule_NK[1]
			nk_k <- suppRule_NK[2]
		}
		
		dims <- lapply(levelObj, function(x) x$dimensions)
		for ( i in 1:length(fullTabObj$Freq) ) {
			cellIndex <- fullTabObj$strID[i]			
			tmpIndex <- which(strObj$strID==cellIndex)	
			# some (sub)totals need to be considered
			if( length(tmpIndex) == 0 ) {
				levInfo <- list()
				for ( z in 1:length(levelObj) ) {
					subLevel <- substr(cellIndex, strObj$strInfo[[z]][1], strObj$strInfo[[z]][2])
					orderInd <- unlist(lapply(dims[[z]], function(x) { match(subLevel, x)}))		
					if( min(orderInd, na.rm=TRUE) == 1 )
						levInfo[[z]] <- dims[[z]][[which(orderInd==1)]]
					else
						levInfo[[z]] <- subLevel	
				}	
				cellIndex <- pasteStrVec(unlist(expand.grid(levInfo)), length(levInfo))
				tmpIndex <- which(strObj$strID %in% cellIndex)	
			}				
			
			contr <- sort(datObj$numVal[tmpIndex], decreasing=TRUE)
			# check for p-percent rule
			if( !is.null(suppRule_P) && length(contr) >= 2 ) 
				suppsP[i] <- f.suppRule.pPerc(sum(contr), contr[1], contr[2], p)						
			# check for nk-rule
			if( !is.null(suppRule_NK) && length(contr) >= nk_n ) 
				suppsNK[i] <- f.suppRule.nk(sum(contr), sum(contr[1:nk_n]), nk_k)	
		}	
	}
	
	indPrimSupps <- which(suppsFreq==TRUE | suppsP==TRUE | suppsNK==TRUE)
	fullTabObj$status[indPrimSupps] <- "u"
	fullTabObj$w[indPrimSupps] <- 0
	outObj$fullTabObj <- fullTabObj
	
	cat("--> Info:", length(which(fullTabObj$status=="u")), "primary cells have been suppressed!\n")
	df <- data.frame(suppRule_Freq=NA, suppRule_P=NA, suppRule_NK=NA)
	df$suppRule_Freq <- sum(suppsFreq)
	df$suppRule_P <- sum(suppsP)
	df$suppRule_NK <- sum(suppsNK)
	rownames(df) <- ""
	cat("--> Info: below you find the distribution of suppressions based on each possible primary suppression rule.\n")
	print(df)
	class(outObj) <- "outObj"
	outObj
}
