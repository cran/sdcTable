# checks if a given vector of secondary suppressions is feasible
checkSuppressionPattern <- function(outObj, pattern, debug=TRUE, all=FALSE, stop=TRUE) {
	# TODO: remove this function and use the one in functions_HITAS_OPT_CHECK!
	attProbOrig <- function(M, subTab, cells2check, xi) {
		A <- M		
		nrVars <- ncol(A)
		nrConstraints <- nrow(A)
		A <- as.simple_triplet_matrix(A)
		dir <- rep("==", nrConstraints)
		rhs <- rep(0, nrConstraints)
		types <- NULL
		
		lb <- subTab$lb
		ub <- subTab$ub
		
		indexSupps <- which(xi!=1)		
		lb[indexSupps] <- subTab$w[indexSupps]
		ub[indexSupps] <- subTab$w[indexSupps]
		
		# we create lpAttProb only once and update it later!	
		ind <- 1:nrVars
		bounds <- list(
				lower=list(ind=ind, lower=lb),
				upper=list(ind=ind, upper=ub)
		)			
		
		limitDown <- limitUp <- rep(NA, length(primSupps))
		for ( i in 1:length(cells2check) ) {
			objF <- rep(0, ncol(A))
			objF[cells2check[i]] <- 1
			
			solUp <- Rglpk_solve_LP(objF, A, dir, rhs, types, max=TRUE, bounds) 
			solDown <- Rglpk_solve_LP(objF, A, dir, rhs, types, max=FALSE, bounds) 
			
			if ( solUp$status == 0 )
				limitUp[i] <- solUp$optimum
			if ( solDown$status == 0 )
				limitDown[i] <- solDown$optimum			
			
		}
		return(list(limitDown=limitDown, limitUp=limitUp))	
	}	
	
	if ( !class(outObj) == "safeTable" )
		stop("'outObj' must be of class 'safeTable'!\n")
	
	if ( !is.vector(pattern) )
		stop("'pattern' must be a vector!\n")	
	
	if ( length(pattern) !=  length(outObj$outObj$strID) )
		stop("please check the length of 'pattern'!\n")
	
	subTab <- as.list(outObj$outObj)
	solver <- "glpk"
	primSupps <- which(subTab$status=="u")
	secondSupps <- which(pattern==TRUE & subTab$status!="u")
	
	xi <- rep(0, length(subTab$strID))
	xi[primSupps] <- 1
	if ( length(secondSupps) > 0 ) {
		xi[secondSupps] <- 1
		subTab$UPL[secondSupps] <- 0
		subTab$LPL[secondSupps] <- 0
		subTab$SPL[secondSupps] <- 0.1
	}
	UPL.p <- subTab$UPL[xi==1]
	LPL.p <- subTab$LPL[xi==1]
	SPL.p <- subTab$SPL[xi==1]
	
	M <- genMatMFull(subTab$strID, outObj$levelObj, check=TRUE)

	# TODO: parallelize
	validPattern <- TRUE
	cells2check <- sort(c(primSupps, secondSupps))
		
	if ( all == FALSE ) {
		primSuppInd <- which(sort(cells2check) %in% primSupps==TRUE)		
		cells2check <- cells2check[primSuppInd]
		UPL.p <- UPL.p[primSuppInd]
		LPL.p <- LPL.p[primSuppInd]
		SPL.p <- SPL.p[primSuppInd]
	}
	
	if ( debug == TRUE )
		cat("Attacker-Problems are now being solved...")
	limits <- attProbOrig(M, subTab, cells2check, xi)
	if ( debug == TRUE )
		cat("[DONE]\n")
	
	limitDown <- limits$limitDown
	limitUp <- limits$limitUp
	
	if ( length(cells2check) > 0 ) {
		if ( all == TRUE & debug == TRUE )
			cat("Checking if all cells (primary and secondary suppressed cells) are adequately protected...\n")
		if ( all == FALSE & debug == TRUE )
			cat("Checking if all primary sensitive cells are adequately protected...\n")
		for ( i in 1:length(cells2check) ) {
			if ( debug == TRUE)
				cat("--> checking cell",i,"|",length(cells2check),"( strID =",subTab$strID[cells2check[i]],"| Freq =",subTab$Freq[cells2check[i]],"): calculated bounds are [",limitDown[i],":",limitUp[i],"]\n")
			checkLimits <- c(UPL.p[i], LPL.p[i], SPL.p[i])
			if ( subTab$Freq[cells2check[i]] + checkLimits[1] > limitUp[i] )
				validPattern <- FALSE
			if ( subTab$Freq[cells2check[i]] - checkLimits[2] < limitDown[i]  )
				validPattern <- FALSE	
			if ( limitUp[i] - limitDown[i] <= checkLimits[3] )
				validPattern <- FALSE
			if ( validPattern == FALSE & stop == TRUE & debug == TRUE ) {
				cat("\nAt least one cell is not safe enough. This pattern is therefore not valid! Stopping here!\n")
				break()
			}				
		}	
		if ( validPattern == TRUE & debug == TRUE )
			cat("\nThis suppression pattern is valid! All cells are adequately protected!\n")
	}
	else {
		if ( debug == TRUE )
			cat("\nThis suppression pattern is valid since no cells have been primary suppressed!\n")
	}
	limits <- list(down=limitDown, up=limitUp)
	return(list(validPattern=validPattern, limits=limits, indices=cells2check))
}
