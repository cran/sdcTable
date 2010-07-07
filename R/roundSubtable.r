roundSubtable <- function(subtab, indexvars=NULL, base=5, maxS=NULL, method=NULL) {
	controlledRounding <- function(subtab, indexvars, base=NULL, maxS=4) {
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
			
			# recode index-variables
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
			mat
		}	
		createBasisR <- function(vals, base) {
			baseR <- rep(base, length(vals))
			baseR[which(vals%%base==0)] <- 0
			baseR[which(is.na(vals))] <- 0
			baseR
		}
		roundDown <- function(vals, base) {
			remaining <- vals %% base
			roundDown <- vals - remaining
			roundDown
		}
		roundUp <- function(rDown, baseR) {
			roundUp <- rDown + baseR
			roundUp
			
		}	
		
		nrVars <- nrow(subtab)
		vals <- subtab$vals
		
		# Recode NA's as 0 to avoid computional problems: quick and dirty!
		NAs <- which(is.na(vals))
		if(length(NAs) > 0) {
			vals[NAs] <- 0
		}
		
		# make sure, all necessary parameters are specified	
		if(is.null(base)) {
			stop("please speficy a base number for controlled rounding procedure!")
		}		
		
		# create the matrix containing addidivity constraints
		M <- createMatM (subtab, indexvars) 
		
		# create basis numbers (if value already a multiple of base, basisR=0!
		basisR <- createBasisR(vals, base)
		
		# round down/up to base
		rDown <- roundDown(vals, base)
		rUp <- roundUp(rDown, basisR) 
		
		# create linear program for controlled rounding procedure
		objF <- con <- dir <- rhs <- xiplus <- ximinus <- NULL
		
		objF <- c(vals, rep(base, 2*nrVars))
		
		for (i in 1:nrow(M)) {
			x <- rep(0, 3*nrVars)	
			
			row <- M[i,]
			indX <- which(row!= 0)		
			indXP 	<- nrVars + indX
			indXM 	<- 2*nrVars + indX
			
			rhs <- c(rhs, -sum(row * rDown))
			rowTimesBase <- row * base
			
			x[indX] <- rowTimesBase[indX]
			x[indXP] <- rowTimesBase[indX]
			x[indXM] <- -rowTimesBase[indX]
			
			con <- rbind(con, x)
			rownames(con) <- NULL
			dir <- c(dir, "==")
		}	
		
		# maximum number of multiples of base-number?
		for (i in 1:nrVars) {
			x <- rep(0, 3*nrVars)	
			x[nrVars+i] <- 1		
			con <- rbind(con, x)
			
			x <- rep(0, 3*nrVars)	
			x[2*nrVars+i] <- 1		
			con <- rbind(con, x)		
			
			rownames(con) <- NULL
			dir <- c(dir, rep("<=",2))
			rhs <- c(rhs, rep(maxS,2))
		}	
		
		# solving the linear problem and returning the rounded table
		erg <- lp("min", objF, con, dir, rhs, binary.vec=1:nrVars, int.vec=(nrVars+1):(3*nrVars))
		solution <- erg$sol
		valsFinal <- NULL
		for (i in 1:nrVars) {
			valsFinal <- c(valsFinal, rDown[i] + base*(solution[i] + solution[nrVars+i] - solution[2*nrVars+i]))
		}
		
		output <- list(valsOrig=vals, valsRounded=valsFinal)
		output
	}
	
	# http://www.blackwell-synergy.com/doi/pdf/10.1111/j.1751-5823.2007.00010.x?cookieSet=1
	randomRounding <- function(subtab, base=NULL) {
		roundUpOrDown <- function(value, base) {
			# value ... cell value
			# floorx ... largest number k with k*base < cell value
			k <- floor(value / base)
			res <- value - k*base
			prob <- res / base
			
			if(value!=0) {			
				# round up or down
				if(runif(1) < prob) {
					value <- (k+1)*base
				}
				else {
					value <- k*base
				}
			}
			value
		}
		
		vals <- subtab$vals
		nrVals <- length(vals)
		
		# remember NA-cells
		NAs <- which(is.na(vals))
		vals[NAs] <- 0
		
		valsFinal <- sapply(vals, roundUpOrDown, base=base)
		
		# NA's should remain NA's
		valsFinal[NAs] <- NA
		
		output <- list(valsOrig=vals, valsRounded=valsFinal)
		output
	}
	
	# simple or rather conventional rounding
	simpleRounding <- function(subtab, base=NULL) {
		roundUpOrDownSimple <- function(value, base) {
			valInBase <- value/base
			reminder <- value %% base
			upOrDown <- valInBase - (value-reminder)/base
			
			if(upOrDown < 0.5) {
				value <- floor(valInBase)*base
			}
			if(upOrDown == 0.5) {
				r <- rbinom(1, 1, prob=0.5)
				ifelse(r==0, value <- floor(valInBase)*base, value <- ceiling(valInBase)*base	)
			}
			if(upOrDown > 0.5) {	
				value <- ceiling(valInBase)*base	
			}
			value
		}
		
		vals <- subtab$vals
		nrVals <- length(vals)
		
		# remember NA-cells
		NAs <- which(is.na(vals))
		vals[NAs] <- 0
	
		valsFinal <- sapply(vals, roundUpOrDownSimple, base=base)
		
		# NA's should remain NA's
		valsFinal[NAs] <- NA
		
		output <- list(valsOrig=vals, valsRounded=valsFinal)
		output		
	}	
	
	### check parameter
	if(is.null(method) | (!method %in% c("controlled", "random", "simple"))) {
		stop("Please speficy a valid method: Choices are simple, random and controlled.")
	}	
	
	# Method: simple Rounding
	if(method=="simple") {
		out <- simpleRounding(subtab, base)
		out <- out$valsRounded
	}
	
	# Method: random Rounding
	if(method=="random") {
		out <- randomRounding(subtab, base)	
		out <- out$valsRounded
	}	
	
	# Method: controlled Rounding
	if(method=="controlled") {
		if(is.null(indexvars)) {
			stop("Please specify the position of the dimensional variables within the given subtable!")
		}
		if(is.null(maxS)) {
			stop("Please specify the maximal number of steps that cell values are allowed to be rounded up or down!")
		}		
		out <- controlledRounding(subtab, indexvars, base, maxS)
		out <- out$valsRounded
	}	
	
	subtab$vals <- out
	subtab
}