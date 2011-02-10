##########################
### preprocessing data ###
##########################
# reduce a csp instance by removing all cells with zero values
reduceSubTable <- function(subTab) {
	idRemoved <- NULL
	restore <- NULL
	if ( any(subTab$Freq == 0) ) {
		idRemoved <- which(subTab$Freq == 0)
		restore <- list()		
		for (i in 1:(length(subTab)-1) ) {
			restore[[i]] <- subTab[[i]][idRemoved]		
			subTab[[i]] <- subTab[[i]][-idRemoved]				
		}
	}		
	return(list(subTab=subTab, idRemoved=idRemoved, restore=restore))
}	

# preprocessing to reduce problem size
# idea: either generate valid constraints for master problem
# or set limits for specific cells already protected to zero
doPreprocessing <- function(xi, aProb, dimM, w, UPL.p, LPL.p, SPL.p, UB, LB, solver, debug=FALSE) {
	primSupps <- which(xi==1)
	checkInd <- rep(TRUE, length(primSupps))		
	resMat <- signVec <- rhsVec <- list()		
	
	for ( i in 1:length(primSupps) ) {
		res <- solve.attProb(aProb, dimM, xi, w, primSupps[i], c(UPL.p[i], LPL.p[i], SPL.p[i]), LB, UB, bridgeless=FALSE, verbose=FALSE, solver)
		if( !is.null(res$constraints$x) ) {
			resMat[[length(resMat)+1]] <- res$constraints[[1]]
			signVec[[length(signVec)+1]] <- res$constraints[[2]]
			rhsVec[[length(rhsVec)+1]] <- res$constraints[[3]]
		}
		else {
			checkInd[i] <- FALSE	
			UPL.p[i] <- LPL.p[i] <- SPL.p[i] <- 0
		}
		if ( debug == TRUE ) {
			cat("checking problematic cell",i,"|",length(primSupps),"... "); flush.console()
			if ( checkInd[i] == TRUE ) 
				cat("[PROBLEM]\n")
			else
				cat("[OK]\n")			
		}		
	}		
	return(list(resMat=resMat, signVec=signVec, rhsVec=rhsVec, UPL.p=UPL.p, LPL.p=LPL.p, SPL.p=SPL.p))
}	
##########################

#########################################################
### generate matrices used in OPT or HITAS procedures ###
#########################################################
# generate Matrix M used in protectHITAS()
genMatM <- function(strID, strInfo) {
	levMat <- matrix(splitStrID(strID, strInfo), nrow=length(strID))
	nrVars <- length(strID)	
	mm <- list()
	
	for( i in 1:length(strInfo) ) {
		f <- apply(levMat[,-i, drop=FALSE], 1, paste, collapse="-")
		spl <- split(1:length(f), f)
		
		if ( length(spl) == nrVars ) {
			mm[[i]] <- matrix(0, nrow=1, ncol=nrVars)
		}
		else {
			mm[[i]] <- matrix(0, nrow=length(spl), ncol=nrVars)
			lenFaktor <- sapply(spl, length)
			for( z in 1:length(spl) ) {
				# slow but working 
				# TODO: check if reordering beforehand is faster
				insert <- rep(1, lenFaktor[z])
				insert[which.min(levMat[spl[[z]],i])] <- -1
				mm[[i]][z, spl[[z]]] <- insert
			}				
		}		
	}	
	mm <- do.call("rbind", mm)
	sums <- apply(mm, 1, function(x) { sum(as.integer(x))})
	if ( any(sums == 0) ) 
		mm <- mm[sums!=0,,drop=F]
	mm
}

# generate Matrix M used in protectOPT()
genMatMFull <- function(strID, levelObj, check=FALSE) {
	if ( check==FALSE )
		m <- do.call("expand.grid", lapply(levelObj, function(x) { x$codesStandard}))
	else
		m <- do.call("expand.grid", lapply(levelObj, function(x) { c(x$codesStandard, x$codeRemoveOrig)}))
	nrVars <- ncol(m)
	nrCells <- nrow(m)
	constraintM <- list()
	f <- apply(m[,,drop=FALSE], 1, function(x) { paste(x, collapse="") })
	order <- match(strID,f)
	m <- m[order,]
	rownames(m) <- 1:nrow(m)
	f <- f[order]		
	
	for ( i in 1:nrVars ) {
		f <- apply(m[,-i, drop=FALSE], 1, function(x) { paste(x, collapse="-") })
		for ( j in 1:length(levelObj[[i]]$dimension) ) {
			x <- subset(m, m[,i] %in% levelObj[[i]]$dimension[[j]])	
			spl <- split(x, f[as.numeric(rownames(x))])
			v <- matrix(0, ncol=nrCells, nrow=length(spl))
			for ( z in 1:length(spl) ) {
				ind <- rep(1,nrow(spl[[z]])) 
				ind[which.min(as.numeric(as.character(spl[[z]][,i])))] <- -1
				v[z, as.integer(rownames(spl[[z]]))] <- ind							
			}
			constraintM[[length(constraintM)+1]] <- v
		}
	}
	constraintM <- do.call("rbind", constraintM)
	return(constraintM)
}		
#########################################################

#################################################################
### work with master-problem for lpSolveAPI and other solvers ###
#################################################################
# create the problem for lpSolveAPI once, update it later
make.master.LPSolve <- function(nrVars, primSupps, forcedCells=NULL, w) {
	lpMaster <- make.lp(0, nrVars)		
	set.objfn(lpMaster, w)
	for ( i in primSupps )
		add.constraint(lpMaster, 1, "=", 1, indices = i)
	
	# force cells for publication
	if ( length(forcedCells) > 0 ) {
		for ( i in forcedCells )
			add.constraint(lpMaster, 1, "=", 0, indices = i)	
	}	
	set.type(lpMaster, 1:length(w), "binary")
	return(lpMaster)		
}	
# create the problem for other solvers once, update it later
make.master.otherSolvers <- function(nrVars, primSupps, forcedCells=NULL, w) {
	objF <- w
	M <- matrix(0, ncol=nrVars, nrow=length(primSupps))
	M[cbind(1:nrow(M),primSupps)] <- 1		
	dir <- rep("==", length(primSupps))
	rhs <- rep(1, length(primSupps))
	
	if ( length(forcedCells) > 0 ) {
		M2 <- matrix(0, ncol=nrVars, nrow=length(forcedCells))	
		M2[cbind(1:length(forcedCells),forcedCells)] <- 1	
		dir2 <- rep("==", length(forcedCells))
		rhs2 <- rep(0, length(forcedCells))	
		M <- rbind(M, M2)
		dir <- c(dir, dir2)
		rhs <- c(rhs, rhs2)
	}	
	types <- rep("B", ncol(M))
	
	M <- as.simple_triplet_matrix(M)
	return(list(objF=objF, M=M, dir=dir, rhs=rhs, types=types))		
}

# update the master problem depending on the solver
update.master <- function(solver, mProb, resMat, signVec, rhsVec) {
	if ( !exists("lpMaster") )
		lpMaster <- NULL		
	
	resMat <- do.call("rbind", resMat)
	signVec <- do.call("c", signVec)
	rhsVec <- do.call("c", rhsVec)
	if ( solver == "lpsolve" ) {
		for ( i in 1:nrow(resMat) ) 
			add.constraint(lpMaster, resMat[i,], signVec[i], rhsVec[i])	
	}
	if ( solver %in% c("glpk", "symphony", "cplex") ) {
		resMat <- as.simple_triplet_matrix(resMat)
		mProb$M <- rbind(mProb$M, resMat)
		mProb$dir <- append(mProb$dir, signVec)
		mProb$rhs <- append(mProb$rhs, rhsVec)
	}	
	if ( solver == "lpsolve" )
		return(lpMaster)
	if ( solver %in% c("glpk", "symphony", "cplex") )
		return(mProb)
}

# solve master problem depending on the solver
solve.master <- function(solver, mProb=NULL) {
	if ( solver != "lpsolve" & is.null(mProb) )
		stop("error in solveMaster!\n")
	if ( ! exists("lpMaster", which(search()=="myGlobalEnv")) ) 
		lpMaster <- NULL
	else 
		lpMaster <- get("lpMaster", pos=which(search()=="myGlobalEnv"))	
	
	if ( !exists("Rcplex") )
		Rcplex <- NULL	
	
	xi <- NULL
	if ( solver == "lpsolve" ) {
		sol <- solve(lpMaster)
		if( sol == 0 ) 
			xi <- get.variables(lpMaster)
	}
	if ( solver == "glpk" ) {
		sol <- Rglpk_solve_LP(
				mProb$objF, 
				mProb$M,
				mProb$dir, 
				mProb$rhs, 
				mProb$types, 
				max = FALSE, 
				NULL, #bounds
				verbose = FALSE)
		if ( sol$status == 0 ) 
			xi <- sol$solution
	}
	if ( solver == "symphony" ) {
		sol <- Rsymphony_solve_LP(
				mProb$objF, 
				mProb$M,
				mProb$dir, 
				mProb$rhs, 
				NULL, #bounds
				mProb$types, 
				max = FALSE)
		if ( sol$status == 0 ) 
			xi <- sol$solution
	}	
	if ( solver == "clpex" ) {
		sense <- rep(NA, mProb$dir)
		sense[mProb$dir=="=="] <- "E"
		sense[mProb$dir=="<="] <- "L"
		sense[mProb$dir==">="] <- "G"
		sol <- Rcplex(
				mProb$objF, 
				mProb$M, 
				mProb$rhs, 
				Qmat = NULL,
				lb = 0, 
				ub = 1, 
				control = list(),
				objsense = "min", 
				sense = sense, 
				vtype = mProb$types, 
				n = 1)	
		if ( sol$status == 0 )  
			xi <- sol$xsol
	}
	return(xi)
}
#################################################################

###################################################################
### work with attacker-problem for lpSolveAPI and other solvers ###
###################################################################
# make attacker-problem for lpSolveAPI
make.attProb.LPSolve <- function(M, freqs) {
	A <- cbind(t(M), diag(length(freqs)), -1*diag(length(freqs)))
	dir <- rep("=", nrow(A))
	lB <- rep("-Inf", nrow(M))
	uB <- rep("Inf", nrow(M))
	
	# we create lpAttProb only once and update it later!
	lpAttP <- make.lp(nrow(A), 0)	
	for( i in 1:ncol(A) ) 
		add.column(lpAttP, A[,i])
	set.constr.type(lpAttP, dir)	
	set.bounds(lpAttP, lower = lB, columns=1:nrow(M))
	set.bounds(lpAttP, upper = uB, columns=1:nrow(M))	
	lpAttP
}

# make attacker-problem for other solvers
make.attProb.otherSolvers <- function(M, freqs) {
	A <- cbind(t(M), diag(length(freqs)), -1*diag(length(freqs)))
	dir <- rep("==", nrow(A))
	rhs <- rep(0, nrow(A))
	lB <- c(rep(-Inf, nrow(M)), rep(0, 2*ncol(M)))
	uB <- c(rep(Inf, nrow(M)), rep(Inf, 2*ncol(M)))	
	objF <- rep(0, ncol(A))
	types <- NULL
	
	lB <- c(rep(-Inf, nrow(M)), rep(0, 2*ncol(M)))
	uB <- c(rep(Inf, nrow(M)), rep(Inf, 2*ncol(M)))			
	# we create lpAttProb only once and update it later!	
	
	A <- as.simple_triplet_matrix(A)
	return(list(objF=objF, A=A, dir=dir, rhs=rhs, types=types, lB=lB, uB=uB))		
}	

# solve the attackers sub-problem (calculate lower|upper bounds) depending on solver
solve.attProb <- function(aProb, dimM, xi, w, primSupp, limits, LB, UB, bridgeless=FALSE, verbose=FALSE, solver) {
	# strengthen constraints
	f.sub.strengthen <- function(para, s) {
		para <- as.numeric(sapply(para, function(x) { min(x, s) } ))
		para
	}	
	
	if ( solver != "lpsolve" & is.null(aProb) )
		stop("error in solve.attProb()!\n")	
	
	if ( ! exists("lpAttP", which(search()=="myGlobalEnv")) ) 
		lpAttP <- NULL
	else 
		lpAttP <- get("lpAttP", pos=which(search()=="myGlobalEnv"))
	
	if ( !exists("Rcplex") )
		Rcplex <- NULL	
	
	c.x <- c.dir <- c.rhs <- NULL
	# imposing upper protection levels
	if ( solver == "lpsolve" ) {
		## prepare input data for dual problem
		objF <- c( rep(0, dimM[1]), w+UB*xi, -(w-LB*xi))
		rhs <- rep(0, length(w))
		rhs[primSupp] <- 1		
	
		set.objfn(lpAttP, objF, 1:length(objF))	
		set.rhs(lpAttP, rhs, 1:length(rhs))	
		
		# imposing upper protection levels
		
		sol <- solve(lpAttP)
		if( sol == 0 ) {
			varsUp <- get.variables(lpAttP)
			objU <- get.objective(lpAttP)
		}
		else 
			stop ("not solvable (up)!\n")	
		
		set.rhs(lpAttP, (-1)*rhs)	
		sol <- solve(lpAttP)
		if( sol == 0 ) {
			varsDown <- get.variables(lpAttP)
			objD <- -get.objective(lpAttP)
		}
		else 
			stop ("not solvable (down)!\n")				
		
	}	
	if ( solver %in% c("glpk", "symphony", "cplex") ) {
		aProb$objF <- c( rep(0, dimM[1]), w+UB*xi, -(w-LB*xi))
		aProb$rhs <- rep(0, length(w))
		aProb$rhs[primSupp] <- 1		
		
		max <- FALSE
		ind <- 1:dim(aProb$A)[1]
		bounds <- list(
				lower=list(ind=ind, val=aProb$lB[ind]),
				upper=list(ind=ind, val=aProb$uB[ind])
		)			
		
		types <- NULL
		if ( solver == "glpk" ) {
			solUp <- Rglpk_solve_LP(aProb$objF, aProb$A, aProb$dir, aProb$rhs, types, max, bounds) 
			solDown <- Rglpk_solve_LP(aProb$objF, aProb$A, aProb$dir, (-1)*aProb$rhs, types, max, bounds) 
			solUp$objval <- solUp$optimum
			solDown$objval <- solDown$optimum
		}
		if ( solver == "symphony" ) {
			solUp <- Rsymphony_solve_LP(aProb$objF, aProb$A, aProb$dir, aProb$rhs, bounds, types, max) 
			solDown <- Rsymphony_solve_LP(aProb$objF, aProb$A, aProb$dir, (-1)*aProb$rhs, bounds, types, max) 
		}
		if ( solver == "cplex" ) {
			sense <- rep(NA, aProb$dir)
			sense[aProb$dir=="=="] <- "E"
			sense[aProb$dir=="<="] <- "L"
			sense[aProb$dir=="<="] <- "G"
			
			solUp <- Rcplex(aProb$objF, aProb$A, aProb$rhs, NULL,aProb$rhs, lb=bounds$lower, ub=bounds$upper, list(), sense, types) 
			solDown <- Rcplex(aProb$objF, aProb$A, (-1)*aProb$rhs, NULL,aProb$rhs, lb=bounds$lower, ub=bounds$upper, list(), sense, types) 
			solUp$solution <- solUp$xopt
			solDown$solution <- solDown$xopt				
			solUp$objval <- solUp$obj
			solDown$objval <- solDown$objval				
		}
		
		if ( solUp$status == 0 ) {
			varsUp <- solUp$solution
			objU <- solUp$objval
		}		
		else 
			stop ("not solvable (up)!\n")				
		if ( solDown$status == 0 ) {
			varsDown <- solDown$solution
			objD <- -solDown$objval										
		}
		else 
			stop ("not solvable (down)!\n")				
	}		
	
	varsUp <- varsUp[(dimM[1]+1):length(varsUp)]
	abUp <- split(varsUp, rep(1:2, each=dimM[2]))	
	cUp <- abUp[[1]]*UB + abUp[[2]]*LB	
	cUp <- f.sub.strengthen(cUp, limits[1])	
	# UPL.p > 0: page 1014 in paper: constraints with rhs == 0 need not to be included!
	if( round(objU,5) < w[primSupp] + limits[1] & bridgeless==FALSE & limits[1] > 0 ) {
		c.x <- rbind(c.x, cUp)
		c.dir <- c(c.dir, ">=")
		c.rhs <- c(c.rhs, limits[1])					
	}			
	
	varsDown <- varsDown[(dimM[1]+1):length(varsDown)]
	abDown <- split(varsDown, rep(1:2, each=dimM[2]))	
	cDown <- abDown[[1]]*UB + abDown[[2]]*LB	
	cDown <- f.sub.strengthen(cDown, limits[2])
	# LPL.p > 0: page 1014 in paper: constraints with rhs == 0 need not to be included!
	if( round(objD,5) > w[primSupp] - limits[2] & bridgeless==FALSE & limits[2] > 0 ) {
		c.x <- rbind(c.x, cDown)
		c.dir <- c(c.dir, ">=")
		c.rhs <- c(c.rhs, limits[2])					
	}				
	
	# imposing sliding protection level
	cSlid <- (abUp[[1]]+abDown[[1]])*UB + (abUp[[2]]+abDown[[2]])*LB	 
	cSlid <- f.sub.strengthen(cSlid, limits[3])	
	# SPL.p > 0: page 1014 in paper: constraints with rhs == 0 need not to be included!
	if( limits[3] > round(objU - objD, 5) & sum(cSlid) > 0 & bridgeless==FALSE & limits[3] > 0 ) {
		c.x <- rbind(c.x, cSlid)
		c.dir <- c(c.dir, ">=")
		c.rhs <- c(c.rhs, limits[3])				
	}				
	if ( bridgeless == TRUE  & round(objU - objD, 5) == 0 ) {
		out <- (abUp[[1]]+abDown[[1]])*UB + (abUp[[2]]+abDown[[2]])*LB
		c.x <- rep(0, length(xi))
		c.x[ out > 0 ] <- 1
		c.x[ primSupp ] <- -1
		c.rhs <- 0
		c.dir <- c(">")		# page 1016 in paper
	}		
	
	if ( verbose == TRUE )
		cat("calculated bounds for cell",primSupp,": [",objD,":",objU,"]\n");flush.console()
	
	add.constraints <- list(x=c.x, dir=c.dir, rhs=c.rhs)
	out <- list(constraints=add.constraints, objD=objD, objU=objU)
	out	
}	

# wrapper to solve attacker problems parallel with snowfall
wrapper.solveAttProb <- function(inObj, aProb, dimM, xi, w, LB, UB, bridgeless, verbose, solver) {
	primSupp <- as.integer(inObj[1])
	limits <- inObj[2:4]
	
	out <- solve.attProb (aProb, dimM, xi, w, primSupp, limits, LB, UB, bridgeless, verbose, solver) 
	return(out)
}	

###################################################################

#############################################################
### create a heuristic solution for the master problem to ###
### speedup overall solution of the master-problem        ###
#############################################################
# set up problem for lpSolve
make.heuristic.lpSolve <- function(M, w, supps, LB, UB) {
	objF <- rep(w, 2)		
	A <- cbind(M, (-1)*M)
	
	dir <- rep("==", nrow(A))
	rhs <- rep(0, nrow(A))		
	types <- rep("C", ncol(A))
	lB <- rep(0, length(UB)+length(LB))
	uB <- c(UB, LB)		
	
	A <- rbind(A, matrix(rep(0, 2*ncol(A)), nrow=2))	
	dir <- c(dir, rep("==",2))				
	rhs <- c(rhs, rep(0,2))	
	
	lpHeuristic <- make.lp(nrow=nrow(A), ncol=length(objF))	
	set.objfn(lpHeuristic, objF, 1:length(objF))
	for( i in 1:ncol(A) ) 
		add.column(lpHeuristic, A[,i])
	set.constr.type(lpHeuristic, dir)	
	set.bounds(lpHeuristic, lower = lB, columns=1:length(objF))
	set.bounds(lpHeuristic, upper = uB, columns=1:length(objF))	
	lpHeuristic
}	

# set up problem for other solvers
make.heuristic.otherSolvers <- function(M, w, supps, LB, UB) {
	objF <- rep(w, 2)		
	A <- cbind(M, (-1)*M)
	
	dir <- rep("==", nrow(A))
	rhs <- rep(0, nrow(A))		
	types <- rep("C", ncol(A))
	lB <- rep(0, length(UB)+length(LB))
	uB <- c(UB, LB)		
	bounds <- list(
			lower=list(ind=1:ncol(A), val=lB),
			upper=list(ind=1:ncol(A), val=uB)
	)		
	
	A <- rbind(A, matrix(rep(0, 2*ncol(A)), nrow=2))	
	dir <- c(dir, rep("==",2))				
	rhs <- c(rhs, rep(0,2))				
	
	A <- as.simple_triplet_matrix(A)
	return(list(objF=objF, A=A, dir=dir, rhs=rhs, types=types, bounds=bounds))		
}

# solve heuristic for lpSolveAPI
solve.heuristic.lpSolve <- function(w, supps, primSupp, limits) {
	if ( !exists("lpHeuristic") )
		lpHeuristic <- NULL

	objF <- rep(w, 2)		
	lenObj <- length(objF)
	
	indSupp1 <- which(supps==1)
	indSupp2 <- indSupp1 + lenObj/2	
	
	set.objfn(lpHeuristic, rep(0, length(indSupp1)), indSupp1)
	set.objfn(lpHeuristic, rep(0, length(indSupp1)), indSupp2)
	i1 <- 1:(lenObj/2)
	i2 <- ((lenObj/2)+1):lenObj
	
	ind <- length(get.rhs(lpHeuristic))
	# upper protection levels
	if ( limits[1] > 0  ) {	
		set.mat(lpHeuristic, ind-1, primSupp, 1)
		set.constr.type(lpHeuristic, "==", ind-1)
		set.rhs(lpHeuristic, limits[1], ind-1)
		
		set.mat(lpHeuristic, ind, lenObj/2 + primSupp, 1)
		set.constr.type(lpHeuristic, "==", ind)
		set.rhs(lpHeuristic, 0, ind)
		sol <- solve(lpHeuristic)	
		if ( sol == 0) {
			solution <- get.variables(lpHeuristic)	
			supps[which(solution[i1]+solution[i2] > 0)] <- 1
		}			
	}		
	
	# lower protection levels
	if ( limits[2] > 0 ) {
		set.mat(lpHeuristic, ind-1, primSupp, 1)
		set.constr.type(lpHeuristic, "==", ind-1)
		set.rhs(lpHeuristic, 0, ind-1)			
		
		set.mat(lpHeuristic, ind, lenObj/2 + primSupp, 1)
		set.constr.type(lpHeuristic, "==", ind)
		set.rhs(lpHeuristic, limits[2], ind)			
		sol <- solve(lpHeuristic)	
		if ( sol == 0) {
			solution <- get.variables(lpHeuristic)	
			supps[which(solution[i1]+solution[i2] > 0)] <- 1
		}				
	}	
	# sliding protection levels
	if ( limits[3] > 0 ) {
		set.mat(lpHeuristic, ind-1, primSupp, 1)
		set.mat(lpHeuristic, ind-1, lenObj/2 + primSupp, 1)
		set.constr.type(lpHeuristic, "==", ind-1)
		set.rhs(lpHeuristic, limits[3], ind-1)					
		
		set.row(lpHeuristic, ind, rep(0, lenObj), indices = 1:lenObj)
		set.constr.type(lpHeuristic, "==", ind)
		set.rhs(lpHeuristic, 0, ind)			
		sol <- solve(lpHeuristic)	
		if ( sol == 0) {
			solution <- get.variables(lpHeuristic)	
			supps[which(solution[i1]+solution[i2] > 0)] <- 1
		}					
	}		
	supps
}

# solve heuristic for other solvers
solve.heuristic.otherSolvers <- function (hProb, supps, primSupp, limits, solver) {
	solveHprob <- function(hProb, supps, solver) { 
		if ( !exists("Rcplex") )
			Rcplex <- NULL			
		if ( solver == "glpk" ) {
			sol <- Rglpk_solve_LP(
					hProb$objF, 
					hProb$A,
					hProb$dir, 
					hProb$rhs, 
					hProb$types, 
					max = FALSE, 
					hProb$bounds, 
					verbose = FALSE)
			if ( sol$status == 0 ) 
				solution <- sol$solution	
		}
		if ( solver == "symphony" ) {
			sol <- Rsymphony_solve_LP(
					hProb$objF, 
					hProb$A,
					hProb$dir, 
					hProb$rhs, 
					hProb$bounds,
					hProb$types, 
					max = FALSE)
			if ( sol$status == 0 ) 
				solution <- sol$solution
		}	
		if ( solver == "clpex" ) {
			sense <- rep(NA, hProb$dir)
			sense[hProb$dir=="=="] <- "E"
			sense[hProb$dir=="<="] <- "L"
			sense[hProb$dir==">="] <- "G"
			sol <- Rcplex(
					hProb$objF, 
					hProb$A, 
					hProb$rhs, 
					Qmat = NULL,
					lb = hProb$bounds$lower$val, 
					ub = hProb$bounds$upper$val, 
					control = list(),
					objsense = "min", 
					sense = sense, 
					vtype = hProb$types, 
					n = 1)	
			if ( sol$status == 0 )  
				solution <- sol$xsol		
		}
		
		i1 <- 1:(ncol(hProb$A)/2)
		i2 <- ((ncol(hProb$A)/2)+1):ncol(hProb$A)
		supps[which(solution[i1]+solution[i2] > 0)] <- 1
		supps				
	}

	indSupp1 <- which(supps==1)
	indSupp2 <- indSupp1 + hProb$A$ncol/2
	hProb$objF[c(indSupp1, indSupp2)] <- 0	
	
	ind <- hProb$A$nrow
	# upper protection levels
	if ( limits[1] > 0  ) {			
		hProb$A <- updateSimpleTripletMatrix(hProb$A, ind-1, primSupp, 1)
		hProb$dir[ind-1] <- "=="
		hProb$rhs[ind-1] <- limits[1]			
		
		hProb$A <- updateSimpleTripletMatrix(hProb$A, ind, length(hProb$objF)/2 + primSupp, 1)
		hProb$dir[ind] <- "=="			
		hProb$rhs[ind] <- 0
		supps <- solveHprob(hProb, supps, solver)		
	}		
	
	# lower protection levels
	if ( limits[2] > 0 ) {
		hProb$A <- updateSimpleTripletMatrix(hProb$A, ind-1, primSupp, 1)
		hProb$dir[ind-1] <- "=="
		hProb$rhs[ind-1] <- 0
		
		hProb$A <- updateSimpleTripletMatrix(hProb$A, ind, length(hProb$objF)/2 + primSupp, 1)
		hProb$dir[ind] <- "=="			
		hProb$rhs[ind] <- limits[2]			
		
		supps <- solveHprob(hProb, supps, solver)	
	}	
	# sliding protection levels
	if ( limits[3] > 0 ) {
		hProb$A <- updateSimpleTripletMatrix(hProb$A, ind-1, primSupp, 1)
		hProb$A <- updateSimpleTripletMatrix(hProb$A, ind-1, length(hProb$objF)/2 + primSupp, -1)
		hProb$dir[ind-1] <- "=="
		hProb$rhs[ind-1] <- limits[3]
		
		hProb$A <- removeFromSimpleTripletMatrix (hProb$A, i=hProb$A$nrow, j=NULL)			
		hProb$dir <- hProb$dir[-ind]
		hProb$rhs <- hProb$rhs[-ind]
		
		supps <- solveHprob(hProb, supps, solver)	
	}			
	supps
}	

# cleanup procedure to remove overprotected cells
cleanup.heuristic <- function(mProb, xi, primSupps, w) {
	if ( class(mProb) == "lpExtPtr" ) {
		nrRows <- length(get.rhs(mProb))
		nrCols <- length(w)
		MM <- matrix(NA, nrow=nrRows, ncol=nrCols)
		for ( i in 1:nrRows ) {
			for ( j in 1:nrCols ) {
				MM[i,j] <- get.mat(mProb, i, j)
			}
		}			
		mDir <- get.constr.type(mProb)
		mRhs <- get.rhs(mProb)
	}
	else {
		MM <- as.matrix(mProb$M)	
		mDir <- mProb$dir
		mRhs <- mProb$rhs
	}
	
	# which suppressions have to be checked
	checkIndices <- setdiff(which(xi==1), primSupps)
	
	if ( length(checkIndices) > 1 ) {
		# order the indices according to weight		
		ordered.checkIndices <- checkIndices[order(w[checkIndices], decreasing=T)]
		delCounter <- 0
		for ( i in 1:length(ordered.checkIndices) ) {
			SUP <- xi
			SUP[ordered.checkIndices[i]] <- 0
			checkMat <- cbind(apply(MM, 1, function(x) { sum(x*SUP)} ), mDir, mRhs)
			result <- apply(checkMat, 1, function(x) { eval(parse(text=paste(x, collapse="")))})
			if ( all(result==TRUE) ) {
				xi <- SUP	
				delCounter <- delCounter +1
			}				
		}				
	}	
	#cat("Total removed redundant suppressions:",delCounter,"!\n")
	xi
}
#############################################################
