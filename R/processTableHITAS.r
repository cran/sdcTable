processTableHITAS <- function(fullData, ub=NULL, lb=NULL, UPLPerc=35, LPLPerc=25, weight="values") {
	# given a subtable, the primary suppressions are protected using the optimal lp formulation (salazaar/fiscetti)
	protectSubtableOpt <- function(subtab, indexvars, lb=NULL, ub=NULL, LPLPerc=NULL, UPLPerc=NULL, weight) {
		# create matrix containing linear constraints
		createMatM <- function(subtab, indexvars) {
			recodeFactor <- function(factorOrig) {
				factorNeu <- rep(0, length(factorOrig))
				origs <- as.character(sort(unique(factorOrig)))
				for (i in 2:length(origs)) {
					factorNeu[which(factorOrig == origs[i])] <- (i-1)
				}
				factorNeu
			}
			
			subtab$indexLP <- 1:nrow(subtab)
			ltab <- length(indexvars)
			nrVars <- nrow(subtab)

			# recode
			for (i in indexvars) {
				subtab[,i] <- recodeFactor(subtab[,i])
			}
			
			anzLevs <- apply(subtab[,indexvars], 2, function(x) { length(unique(x)) } )
			mat <- NULL
			for (z in 1:length(indexvars)) {
				# create factor from other indexvars
				ind <- indexvars[-z]
				cmd <- "factor(paste("
				for (j in 1:length(ind)) {
					cmd <- paste(cmd, "subtab[,", ind[j],"], ", sep="")
				}
				cmd <- paste(cmd, "sep=\"\"))")
				subtab$fac <- eval(parse(text=cmd))

				# split data and calculate matrices and vectors for linear optimization programm
				t <- split(subtab, subtab$fac)
				for(i in 1:length(t)) {
					x <- rep(0, nrVars)
					if(t[[i]][1,z] == min(t[[i]][,z])) {
						x[t[[i]]$indexLP] <- c(-1, rep(1, nrow(t[[i]])-1))
					}
					else {
						x[t[[i]]$indexLP] <- c(rep(1, nrow(t[[i]])-1), -1)
					}
					mat <- rbind(mat, x)
				}
			}
			rownames(mat) <- NULL
			mat <- mat[sample(1:nrow(mat)),]
		}
		# the Attacker Subproblem (minimize, maximize)
		AttackerSubproblem <- function(vals, xi, M, primSupp, LB, UB) {
			
			# quick and dirty fix for NA-values
			vals[which(is.na(vals))] <- 1
			
			con <- M
			con <- rbind(M, -1*M)
			dir <- rep("<=", nrow(con))
			rhs <- rep(0, nrow(con))

			# y_i <= a_i + UB_i
			con <- rbind(con, diag(ncol(M)))
			dir <- c(dir, rep("<=", ncol(M)))
			rhs <- c(rhs, vals+UB*xi)

			# -y_i <= -(a_i - LB_i)
			con <- rbind(con, -1*diag(ncol(M)))
			dir <- c(dir, rep("<=", ncol(M)))
			rhs <- c(rhs, -(vals-LB*xi))

			objF <- rep(0, length(vals))
			objF[primSupp] <- 1

			# LP- Duality
			rhsDual <- objF
			objFDual <- rhs
			conDual <- t(con)
			dirDual <- rep(">=", length(rhsDual))

			anzVars <- nrow(diag(ncol(conDual)))

			# gamma_i, alpha_i, beta_i >= 0
			conDual <- rbind(conDual, diag(ncol(conDual)))
			dirDual <- c(dirDual, rep(">=", anzVars))
			rhsDual <- c(rhsDual, rep(0, anzVars))

			ergDualU <- lp("min", objFDual, conDual, dirDual, rhsDual)
			ergDualD <- lp("min", objFDual, conDual, dirDual, -1*rhsDual)
			limitDown <- -ergDualD$objv
			limitUp <- ergDualU$objv
			return(list(ergDown=ergDualD, ergUp=ergDualU, limitDown=limitDown, limitUp=limitUp, vals=vals))
		}
		# solving the master Linear Program
		solveMasterLP <- function(vals, primSupps, res=NULL) {
			anzPrimSupps <- length(primSupps)
			anzVars <- length(vals)

			# objective Function
			objF <- vals

			# constraints
			con <- matrix(0, ncol=anzVars, nrow=anzPrimSupps)
			dir <- NULL
			rhs <- NULL
			for (i in 1:anzPrimSupps) {
				con[i, primSupps[i]] <- 1
				dir <- c(dir, "==")
				rhs <- c(rhs, 1)
			}

			# add additional restrictions (from solutions in attackers subproblems)
			if(!is.null(res)) {
				con <- rbind(con, res$con)
				rhs <- c(rhs, res$rhs)
				dir <- c(dir, res$dir)
			}

			erg <- lp("min", objF, con, dir, rhs, binary.vec=1:length(objF))
		    erg
		}
		# combine and (possibly) add new restrictions
		combineRestrictions <- function(resOld, attackerProblem, LB, UB, LPL, UPL, primSupp) {
			strenghtenConstraints <- function(con, rhs) {
				con <- as.numeric(sapply(con, function(x) { min(x, rhs) } ))
				con
			}
			coverInequalities <- function(con) {
				conInequality <- rep(0, length(conNew))
				conInequality[which(con != 0)] <- 1
				conInequality
			}

			solLow 	<- attackerProblem$ergDown$sol
			solUp 	<- attackerProblem$ergUp$sol
			boundLow <- attackerProblem$limitDown
			boundUp <- attackerProblem$limitUp
			vals 	<- attackerProblem$vals
			anzVars <- length(solLow)	

			endBeta 	<- anzVars
			startBeta	<- anzVars - length(UPL) + 1
			endAlpha 	<- startBeta - 1
			startAlpha 	<- endAlpha - length(UPL) + 1	

			ind <- 0	
			con <- rhs <- dir <- NULL

			if(boundUp < vals[primSupp] + UPL[primSupp]) {
				conNew <- (solUp[startAlpha:endAlpha] * UB) + (solUp[startBeta:endBeta]*LB)
				rhsNew <- UPL[primSupp]
				dirNew <- ">="

				# Strenghten Capacity constraints
				conNew <- strenghtenConstraints(conNew, rhsNew)
				con <- rbind(con, conNew)
				dir <- c(dir, dirNew)
				rhs <- c(rhs, rhsNew)

				# cover inequalities
				conInequality <- coverInequalities(conNew)
				con <- rbind(con, conInequality)
				dir <- c(dir, ">=")
				rhs <- c(rhs, 1) 
			}

			if(boundLow > vals[primSupp] - LPL[primSupp]) {
				conNew <- (solLow[startAlpha:endAlpha] * UB) + (solLow[startBeta:endBeta]*LB)
				rhsNew <- LPL[primSupp]
				dirNew <- ">="

				# Strenghten Capacity constraints
				conNew <- strenghtenConstraints(conNew, rhsNew)

				con <- rbind(con, conNew)
				dir <- c(dir, dirNew)
				rhs <- c(rhs, rhsNew)

				# cover inequalities
				conInequality <- coverInequalities(conNew)
				con <- rbind(con, conInequality)
				dir <- c(dir, ">=")
				rhs <- c(rhs, 1) 
			}
			
			con <- rbind(resOld$con, con)
			
			dir <- c(resOld$dir, dir)
			rhs <- c(resOld$rhs, rhs)
			
			res <- list(con=con, dir=dir, rhs=rhs)	
			res
		}

		nrVars <- nrow(subtab)		
		fixedVals <- subtab$val
		fixedVals[is.na(fixedVals)] <- 1		
		
		# lb = lower bound
		if(is.null(lb)) {
			lb <- rep(0, nrVars)
		}

		# lb = lower bound
		if(is.null(ub)) {
			ub <- fixedVals * 100
		}

		UPL <- ceiling(fixedVals * (1+(UPLPerc/100))) - fixedVals 
		UPL <- as.numeric(sapply(UPL, function(x) { max(x, 1) } ))

		LPL <- fixedVals - ceiling(fixedVals * ((100-LPLPerc)/100)) 
		LPL <- as.numeric(sapply(LPL, function(x) { max(x, 1) } ))

		UB <- ub - fixedVals
		LB <- fixedVals - lb

		primSupp <- which(subtab$geh == "P")
		indNAS <- which(is.na(subtab$val))
		if(length(indNAS) > 0) {
			primSupp <- c(primSupp, indNAS)
		}
		
		# other Suppressions
		secondarySupps <- which(subtab$geh != "" & subtab$geh != "P")
		
		xi <- rep(0, nrow(subtab))
		xi[primSupp] <- 1
		allSuppsStart <- primSupp
		
		if(length(secondarySupps) > 0) {
			xi[secondarySupps] <- 1
			allSuppsStart <- c(allSuppsStart, secondarySupps)
		}

		# create the matrix containing addidivity constraints
		M <- createMatM (subtab, indexvars) 

		ind <- FALSE
		res <- NULL
		nrConstraints <- 0

		# use weight for objective-function
		# values: w_i equals cell values
		# logs: w_i equals log(1+cell values)
		if(!weight %in% c("values","logs")) {
			stop("You need to speficy a correct weighting scheme. Possible choices are: values, logs.")
		}
		else {
			indSupps <- which(subtab$geh %in% c("S", "P"))
			notindSupps <- which(subtab$geh =="")
			v <- rep(NA, nrVars)
			v[indSupps] <- 0
			if(weight=="values") {					
				v[notindSupps] <- fixedVals[notindSupps]	
			}
			if(weight=="logs") {
				v[notindSupps] <- log(1+fixedVals[notindSupps])
			}					
		}			
		
		xi <- solveMasterLP(v, allSuppsStart, res=NULL)$sol
		while (ind == FALSE) {	
			if(!is.null(res)) { 	
				nrConstraints <- nrow(res$con)
			}
			for (z in 1:length(primSupp)) {
				attProb <- AttackerSubproblem(subtab$val, xi, M, primSupp[z], LB, UB)
				res <- combineRestrictions(res, attProb, LB, UB, LPL, UPL, primSupp[z])
			}
			
			xi <- solveMasterLP(v, allSuppsStart, res)$sol
			if(nrConstraints == nrow(res$con)) { 	
				ind <- TRUE
			}	
		}

		# recode subtab: S... secondary suppressed cells, P ... primary suppressed cells
		subtab$geh[which(xi==1)] <- "S"
		subtab$geh[primSupp] <- "P"
		secondarySupps <- as.numeric(rownames(subtab)[which(subtab$geh=="S")])
		secondarySupps
	}

	start <- as.numeric(Sys.time())
	
	# primary Supps are "fullData$supps2check", all other supps are "secondary"
	SuppsStart <- which(fullData$data$geh != "")	
	if(length(SuppsStart) > 0) {
		fullData$data$geh[SuppsStart] <- "S"
	}
	fullData$data$geh[fullData$supps2check] <- "P"
	
	anzSuppStart <- length(SuppsStart)
	spl <- splitPrimarySupps(fullData)
	nc <- ncol(spl[[1]])

	cat("The algorithm runs over all subgroups with primary/secondary suppressed values!\n")
	lSpl <- length(spl)
	for (i in 1:lSpl) {
		txtProgressBar(max = lSpl, initial = i, char = "=",style = 3)
		subtab <- calcSubset(fullData, spl[[i]])

		suppsInFullData <- protectSubtableOpt(subtab$data, fullData$indexvars, lb, ub, LPLPerc, UPLPerc, weight)
		if(length(suppsInFullData) > 0) {
			fullData$data[suppsInFullData, "geh"] <- "S"			
			# is it only necessary to split the dataset again if new suppressions are added to the data
			spl <- splitPrimarySupps(fullData)
		}
	}
	cat("\n")

	SuppsEnd <- which(fullData$data$geh!= "")
	anzSuppEnd <- length(SuppsEnd)

	# all cells with NA's are also secondary suppressions
	indNA <- which(is.na(fullData$data$val))
	if(length(indNA) > 0) {
		fullData$data$geh[indNA] <- "S"
	}
	rm(indNA)	
	
	suppsNew <- SuppsEnd[-which(SuppsEnd %in% SuppsStart)]

	anzSuppSec <- anzSuppEnd - anzSuppStart
	end <- as.numeric(Sys.time())

	fullData$supps2check <- suppsNew

    erg <- list()
    erg$fullData <- fullData
	erg$counter <- 1
    erg$time <- end-start
	erg$method <- "OPT"
	erg$totSupps <- length(which(fullData$data$geh=="S"))
	
	class(erg) <- "safeTable"
    return(erg)

	#return(list(fullData=fullData, anzSecSupp=anzSuppSec, time=end-start))
}