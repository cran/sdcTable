protectLinkedTables <- function(inputObj1, inputObj2, commonCells, method="HITAS", weight=NULL) {
	f.calcCommonCellIndices <- function(outObj1, levelObj1, outObj2, levelObj2, commonCells) {
		commonInd1 <- 1:length(outObj1$strID)
		commonInd2 <- 1:length(outObj2$strID)
		
		varsUsed1 <- as.numeric(unique(unlist(lapply(commonCells, function(x) x[1]))))
		varsUsed2 <- as.numeric(unique(unlist(lapply(commonCells, function(x) x[2]))))
		varsNotUsed1 <- setdiff(1:length(levelObj1), varsUsed1)
		varsNotUsed2 <- setdiff(1:length(levelObj2), varsUsed2)
		
		startVar <- which(names(res1)=="LB")
		
		for ( i in 1:length(commonCells) ) {
			# it is not te same variable --> different characterisics
			if ( length(commonCells[[i]]) != 3 ) {
				commonInd1 <- setdiff(commonInd1, which(!outObj1[[startVar+as.numeric(commonCells[[i]][[1]])]] %in% commonCells[[i]][[3]]))			
				commonInd2 <- setdiff(commonInd2, which(!outObj2[[startVar+as.numeric(commonCells[[i]][[2]])]] %in% commonCells[[i]][[4]]))			
			}		
		}	
		
		if ( length(varsNotUsed1) > 0 ) {
			for ( i in varsNotUsed1 ) 			
				commonInd1 <- setdiff(commonInd1, which(outObj1[[startVar+i]] != levelObj1[[c(varsNotUsed1[i])]]$codesStandard[1]))				
		}
		
		if ( length(varsNotUsed2) > 0 ) {
			for ( i in varsNotUsed2 ) 			
				commonInd2 <- setdiff(commonInd2, which(outObj2[[startVar+i]] != levelObj2[[c(i)]]$codesOrig[1]))				
		}
		
		if ( length(commonInd1) != length(commonInd2) )
			stop("generic error in calcCommonCellIndices!\n")
		return(list(commonInd1=commonInd1, commonInd2=commonInd2))
	}				
	
	f.checkCommonCells <- function(suppPattern1, suppPattern2, commonCellIndices) {
		indOK <- TRUE
		if ( any(suppPattern1[commonCellIndices[[1]]] != suppPattern2[commonCellIndices[[2]]]) ) 
			indOK <- FALSE
		return(indOK)
	}
	
	if ( !method %in% c("HYPERCUBE","HITAS") )
		stop("please specify a suitable protection method!\n")
		
	levelObj1 <- inputObj1$levelObj
	levelObj2 <- inputObj2$levelObj
	
	# firstRun
	outObj1 <- protectTable(outObj=inputObj1, method=method)		
	outObj2 <- protectTable(outObj=inputObj2, method=method)		
	
	res1 <- outObj1$outObj
	res2 <- outObj2$outObj
	
	# calc original primary suppressions
	origPrimSupp1Index <- which(res1$status == "u")
	origPrimSupp2Index <- which(res2$status == "u")
	
	# calculate commonCells:
	commonCellIndices <- f.calcCommonCellIndices(res1, levelObj1, res2, levelObj2, commonCells)	
	
	suppPattern1 <- rep(0, length(res1$status))
	suppPattern1[which(res1$status %in% c("u","x"))] <- 1
	suppPattern2 <- rep(0, length(res2$status))
	suppPattern2[which(res2$status %in% c("u","x"))] <- 1	
	
	indOK <- f.checkCommonCells(suppPattern1, suppPattern2, commonCellIndices)
	if ( indOK == FALSE ) {
		runInd <- TRUE
		counter <- 1
		while ( runInd == TRUE ) {
			x <- cbind(suppPattern1[commonCellIndices[[1]]], suppPattern2[commonCellIndices[[2]]])
			index <- list()
			i1 <- which(x[,1] == 0 & x[,2]==1)
			i2 <- which(x[,1] == 1 & x[,1]==0)
			index[[1]] <- commonCellIndices[[1]][i1]
			index[[2]] <- commonCellIndices[[2]][i2]
			
			for ( j in 1:2 ) {
				if ( length(index[[j]]) > 0 ) {
					if ( j == 1 ) {
						inputObj1$fullTabObj$status[na.omit(match(outObj1$outObj$strID[index[[j]]], inputObj1$fullTabObj$strID))] <- "u"
						outObj1 <- protectTable(outObj=inputObj1, method=method)	
					}
					else {
						inputObj2$fullTabObj$status[na.omit(match(outObj2$outObj$strID[index[[j]]], inputObj2$fullTabObj$strID))] <- "u"
						outObj2 <- protectTable(outObj=inputObj2, method=method)	
					}					
				}				
			}
			suppPattern1 <- rep(0, length(outObj1$outObj$status))
			suppPattern1[which(outObj1$outObj$status %in% c("u","x"))] <- 1
			suppPattern2 <- rep(0, length(outObj2$outObj$status))
			suppPattern2[which(outObj2$outObj$status %in% c("u","x"))] <- 1				
						
			cbind(suppPattern1[commonCellIndices[[1]]], suppPattern2[commonCellIndices[[2]]])
			indOK <- f.checkCommonCells(suppPattern1, suppPattern2, commonCellIndices)
			if ( indOK == TRUE)
				runInd <- FALSE	
			if( counter > 10) {
				runInd <- FALSE
				stop("iterative procedure did not converge.\n")
			}				
			counter <- counter + 1			
		}		
	} else { 
		cat("\n===> all common cells have the same anonymity-state in both tables! [Finished]\n")
	}
	# return objects	
	#ind1 <- which(outObj1$outObj$primSupp[-origPrimSupp1Index]==TRUE)
	#ind2 <- which(outObj2$outObj$primSupp[-origPrimSupp2Index]==TRUE)
	#if ( length(ind1) > 0 )
	#	outObj1$outObj$secondSupps[ind1] <- FALSE
	#if ( length(ind2) > 0 )
	#	outObj2$outObj$secondSupps[ind2] <- FALSE	
	return(list(outObj1=outObj1, outObj2=outObj2))
}
