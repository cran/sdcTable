# TODO: implement parameter "weight" for HITAS procedure
# currently, FREQ is used by default!
protectHITAS <- function(subTab, M, solver, verbose=FALSE, cpus=1, method) {
	error <- FALSE
	Rcplex <- Rsymphony <- Rglpk <- NULL	
	lpHeuristic <- lpAttP <- lpMaster <- NULL
	sfInit <- sfLibrary <- sfExport <- sfClusterApplyLB <- sfStop <- NULL
	#strInfo <- strObj$strInfo
	
	#tt <- proc.time()
	#cat("### reducing the subTable...")
	#idRemoved <- NULL
	#if ( method != "OPT" ) {
	#	red <- reduceSubTable(subTab)
	#	subTab <- red$subTab
	#	idRemoved <- red$idRemoved
	#	restore <- red$restore
	#	rm(red)			
	#}	
	#cat(" [DONE in",formatC(round((proc.time()-tt)[3]/60,2), digits=2, format="f"),"minutes] ###\n")
	
	primSupps <- which(subTab$status=="u")
	secSuppsInd <- which(subTab$status=="x")
	forcedCells <- which(subTab$status=="z")
	UPL.p <- LPL.p <- SPL.p <- NULL
	
	UPL.p <- subTab$UPL[primSupps]
	LPL.p <- subTab$LPL[primSupps]
	SPL.p <- subTab$SPL[primSupps]
	
	if ( length(secSuppsInd) > 0 ) {
		primSupps <- c(primSupps, secSuppsInd)
		UPL.p <- c(UPL.p, rep(0, length(secSuppsInd)))
		LPL.p <- c(LPL.p, rep(0, length(secSuppsInd)))
		SPL.p <- c(SPL.p, rep(1, length(secSuppsInd)))
	}
	
	dimM <- dim(M)
	
	# set up attacker and master problems
	#tt <- proc.time()
	#cat("### setting up the attackers and master problem...")	
	mProb <- aProb <- NULL
	lpAttP <- lpMaster <- NULL
	if ( solver %in% c("glpk", "symphony", "cplex") ) {
		mProb <- make.master.otherSolvers(ncol(M), primSupps, forcedCells, subTab$w)	
		aProb <- make.attProb.otherSolvers(M, subTab$Freq)
	}
	if ( solver == "lpsolve" ) {
		lpMaster <- make.master.LPSolve(ncol(M), primSupps, forcedCells, subTab$w)	
		lpAttP <- make.attProb.LPSolve(M, subTab$Freq)
		assign("lpMaster", lpMaster, pos=which(search()=="myGlobalEnv"))
		assign("lpAttP", lpAttP, pos=which(search()=="myGlobalEnv"))
	}		
	#cat(" [DONE in",round((proc.time()-tt)[3]/60,2),"minutes] ###\n")
		
	# we start by suppressing all primary and secondary suppressed cells
	xi <- rep(0, length(subTab$w))
	xi[c(primSupps, secSuppsInd)] <- 1
	
	#if ( method == "OPT" ) {
	#	tt <- proc.time()
	#	#cat("### start pre-processing...")	
	#	res.pre <- doPreprocessing(xi, aProb, dimM, subTab$w, UPL.p, LPL.p, SPL.p, subTab$UB, subTab$LB, solver, debug=FALSE) 
	#	UPL.p <- res.pre$UPL.p
	#	LPL.p <- res.pre$LPL.p
	#	SPL.p <- res.pre$SPL.p
	#	
	#	if ( !length(res.pre$resMat)==0 )
	#		mProb <- update.master(solver, mProb, res.pre$resMat, res.pre$signVec, res.pre$rhsVec) 
	#	rm(res.pre)
	#	#cat(" [DONE in",formatC(round((proc.time()-tt)[3]/60,2), digits=2, format="f"),"minutes] ###\n")
	#	
	#	### generate heuristic solution
	#	#cat("### generating a heuristic solution...")	
	#	if ( solver == "lpsolve" ) {
	#		cat("### generate heuristic lp-Problem with lpSolve...")
	#		lpHeuristic <- make.heuristic.lpSolve(M, subTab$w, xi, subTab$LB, subTab$UB)
	#		cat("[DONE]\n")
	#		for ( i in 1:length(primSupps) ) {
	#			limits <- c(UPL.p[i], LPL.p[i], SPL.p[i])
	#			xi <- solve.heuristic.lpSolve(subTab$w, xi, primSupps[i], limits)				
	#			#cat("nrSupps in run",i,":", sum(xi),"\n")
	#		}			
	#	}
	#	if ( solver %in% c("glpk","symphony","cplex") ) {
	#		hProb <- make.heuristic.otherSolvers(M, subTab$w, xi, subTab$LB, subTab$UB)	
	#		for ( i in 1:length(primSupps) ) {
	#			limits <- c(UPL.p[i], LPL.p[i], SPL.p[i])
	#			xi <- solve.heuristic.otherSolvers(hProb, xi, primSupps[i], limits, solver)
	#			#cat("nrSupps in run",i,":", sum(xi),"\n")
	#		}				
	#	}	
	#	# do cleanup
	#	xi <- cleanup.heuristic(mProb, xi, primSupps, subTab$w)
	#	#cat(" [DONE in",formatC(round((proc.time()-tt)[3]/60,2), digits=2, format="f"),"minutes] ###\n")
	#}
	
	runInd <- FALSE
	counter <- 1	
	while ( runInd == FALSE ) {
		#if ( counter %% 10 == 1 )
		#	cat("we are starting Run",counter,"\n"); flush.console()
		#cat("we are now solving attacker-subproblems!\n")
		resMat <- signVec <- rhsVec <- list()		
		limitUp <- limitDown <- rep(NA, length(primSupps))
		#tt <- proc.time()
		#cat("### solving attacker problems...")			
	
		# we can't do parallelization
		if ( cpus==1 | solver=="lpsolve" ) {
			for ( i in 1:length(primSupps) ) {
				limits <- c(UPL.p[i], LPL.p[i], SPL.p[i])
				if ( any(limits > 0) ) {
					res <- solve.attProb(aProb, dimM, xi, subTab$Freq, primSupps[i], limits, subTab$LB, subTab$UB, bridgeless=FALSE, verbose=FALSE, solver)
					
					out <- res$constraints
					limitDown[i] <- res$objD
					limitUp[i] <- res$objU
					#cat("  --> Cell",i,"( ID:",subTab$strID[primSupps[i]],"): [",limitDown[i],":",limitUp[i],"]\n"); flush.console()
					if( !is.null(out$x) ) {
						resMat[[length(resMat)+1]] <- out[[1]]
						signVec[[length(signVec)+1]] <- out[[2]]
						rhsVec[[length(rhsVec)+1]] <- out[[3]]
					}
				}						
			}			
		}
		else {
			#cat("\n==> Note: starting parallel processing with",cpus,"processors!\n")
			inObj <- list()
			for ( i in 1:length(primSupps) ) {
				limits <- c(UPL.p[i], LPL.p[i], SPL.p[i])
				if ( any(limits > 0) )
					inObj[[length(inObj)+1]] <- c(primSupps[i], limits)	
			}						
			sfInit(parallel=TRUE, cpus=cpus)
			if ( solver == "glpk" )
				sfLibrary(Rglpk, verbose=FALSE)
			if ( solver == "symphony" )
				sfLibrary(Rsymphony, verbose=FALSE)
			if ( solver == "cplex" )
				sfLibrary(Rcplex)
			sfExport("solve.attProb")
			sfExport("wrapper.solveAttProb")			
			res <- sfClusterApplyLB(inObj, wrapper.solveAttProb, aProb, dimM, xi, subTab$Freq, subTab$LB, subTab$UB, bridgeless=FALSE, verbose=FALSE, solver)
			sfStop()
			resMat <- lapply(res, function(x) { x$constraints$x })
			signVec <- lapply(res, function(x) { x$constraints$dir })
			rhsVec <- lapply(res, function(x) { x$constraints$rhs })
			limitDown <- sapply(res, function(x) { x$objD })
			limitUp <- sapply(res, function(x) { x$objU })
			
			if (counter == 1) 
				checkInd[sapply(resMat, is.null)] <- FALSE	
		}		
		#cat(" [DONE in",formatC(round((proc.time()-tt)[3]/60,2), digits=2, format="f"),"minutes] ###\n")
				
		#tt <- proc.time()
		#cat("### searching for bridgeless inequalities...")		
		candidates <- setdiff(which(xi == 1), primSupps)
		if( length(candidates) > 0 ) {
			if ( cpus==1 | solver=="lpsolve" ) {
				for (i in 1:length(candidates)) {
					res <- solve.attProb(aProb, dimM, xi, subTab$Freq, candidates[i], c(0, 0, 0.1), subTab$LB, subTab$UB, bridgeless=TRUE, verbose=FALSE, solver)
					br <- res$constraints
					if ( !is.null(br[[1]]) ) {
						resMat[[length(resMat)+1]] <- br[[1]]
						signVec[[length(signVec)+1]] <- br[[2]]
						rhsVec[[length(rhsVec)+1]] <- br[[3]]					
					}
				}					
			}
			else {
				inObj <- list()
				for ( i in 1:length(candidates) ) 
					inObj[[length(inObj)+1]] <- c(candidates[i], c(0, 0, 0.1))	
										
				sfInit(parallel=TRUE, cpus=cpus)
				if ( solver == "glpk" )
					sfLibrary(Rglpk, verbose=FALSE)
				if ( solver == "symphony" )
					sfLibrary(Rsymphony, verbose=FALSE)
				if ( solver == "cplex" )
					sfLibrary(Rcplex)
				sfExport("solve.attProb")
				sfExport("wrapper.solveAttProb")			
				res <- sfClusterApplyLB(inObj, wrapper.solveAttProb, aProb, dimM, xi, subTab$Freq, subTab$LB, subTab$UB, bridgeless=TRUE, verbose=FALSE, solver)
				sfStop()
				brMat <- lapply(res, function(x) { x$constraints$x })
				brsignVec <- lapply(res, function(x) { x$constraints$dir })
				brrhsVec <- lapply(res, function(x) { x$constraints$rhs })
				if ( any(sapply(brsignVec, function(x) { is.null(x)})==FALSE) )	 {
					resMat <- append(resMat, brMat)
					signVec <- append(signVec, brsignVec)
					rhsVec <- append(rhsVec, brrhsVec)
				}
			}	
		}
		#cat(" [DONE in",formatC(round((proc.time()-tt)[3]/60,2), digits=2, format="f"),"minutes] ###\n")
				
		if ( length(which(sapply(resMat, function(x) {!is.null(x) } ) == TRUE )) > 0 ) {
			#tt <- proc.time()
			#cat("### adding constraints to the master LP...")
			mProb <- update.master(solver, mProb, resMat, signVec, rhsVec)
			#cat(" [DONE in",formatC(round((proc.time()-tt)[3]/60,2), digits=2, format="f"),"minutes] ###\n")
			
			#tt <- proc.time()
			#cat("### solving the master LP...")			
			xi <- solve.master(solver, mProb)
			#cat(" [DONE in",formatC(round((proc.time()-tt)[3]/60,2), digits=2, format="f"),"minutes] ###\n")
			
			if ( is.null(xi) ) {
				#warning("Master-LP cannot be solved!\n")	
				runInd <- TRUE	
				subTab <- NULL
				error <- TRUE
			}
				
			#cat("====> xi in run", counter, "has",sum(xi),"suppressions\n"); flush.console()
			counter <- counter + 1
		}
		else 
			runInd <- TRUE	
	}
	#cat("\n==> calculated bounds for primary suppressed cells in the current (sub)table:\n"); flush.console()
	
	if ( any(is.na(limitDown)) & error == FALSE ) {
		cat("problem with NA's in limitDown...\n")
		for ( i in 1:length(primSupps) ) {
			limits <- c(subTab$UPL[primSupps[i]], subTab$LPL[primSupps[i]], subTab$SPL[primSupps[i]])
			res <- solve.attProb(aProb, dimM, xi, subTab$Freq, primSupps[i], limits, subTab$LB, subTab$UB, bridgeless=FALSE, verbose=FALSE, solver)
			limitDown[i] <- res$objD
			limitUp[i] <- res$objU
		}					
	}
	#for ( i in 1:length(primSupps) )
	#	cat("  --> Cell",i,"( ID:",subTab$strID[primSupps[i]],"): [",limitDown[i],":",limitUp[i],"]\n")

	if ( error == FALSE ) {
		suppIndex <- setdiff(which(xi==1), which(subTab$status=="u"))
		#subTab$secondSupps[suppIndex] <- TRUE
		subTab$status[suppIndex] <- "x"		
	}
	return(subTab)		
}
