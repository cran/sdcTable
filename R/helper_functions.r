## Wrapper function for pasting key-Variables
pasteStrVec <- function(strVec, nrVars, coll=NULL) {
	if(length(strVec) %% nrVars != 0)
		stop("Wrong Dimensions!\n")
	
	if ( is.null(coll) ) {
		.Call( "myPaste", as.character(strVec), nrVars, PACKAGE = "sdcTable")
	} else {
		.Call( "myPasteWithSep", as.character(strVec), nrVars, coll, PACKAGE = "sdcTable")
	}
}

# alternative to expand.grid (used for pasteStrVec!)
expand <- function(inputList, vector=TRUE) {
	uniques <- sapply(inputList, length)
	nrPoss <- prod(uniques)
	if ( vector == TRUE ) {
		out <- NULL
		for ( i in 1:length(inputList) ) {
			if ( i == 1 ) 
				out <- rep(inputList[[i]], nrPoss/length(inputList[[i]]))
			else 
				out <- c(out, rep(inputList[[i]], each=prod(uniques[1:(i-1)]), nrPoss/length(rep(inputList[[i]], each=prod(uniques[1:(i-1)])))))
		}		
	}
	else {
		out <- list()
		for ( i in 1:length(inputList) ) {
			if ( i == 1 ) 
				out[[i]] <- rep(inputList[[i]], nrPoss/length(inputList[[i]]))
			else 
				out[[i]] <- rep(inputList[[i]], each=prod(uniques[1:(i-1)]), nrPoss/length(rep(inputList[[i]], each=prod(uniques[1:(i-1)]))))
		}		
	}
	out
}	

# returns a vector original size or str
mySplit <- function(strVec, keepIndices) {
	if ( min(keepIndices) < 1 | max(keepIndices) > nchar(strVec[1]) ) {
		stop("indices must be in 1:",nchar(strVec[1]),"!\n")	
	}	
	keepIndices <- unique(keepIndices)-1 # required because of indexing in c++
	return(.Call( "mySplitFn", as.character(strVec), as.numeric(keepIndices), PACKAGE = "sdcTable"))
}

#strs <- rep(paste(LETTERS[1:6],collapse=""), 10000)
#system.time({
#	sapply(strs, mySplit, c(1,6))			
#})

mySplitIndicesList <- function(strVec, keepList, coll="-") {
	u <- unlist(keepList)
	if ( min(u) < 1 | max(u) > nchar(strVec[1]) ) {
		stop("indices must be in 1:",nchar(strVec[1]),"!\n")	
	}		
	out <- list()
	for ( i in 1:length(keepList) ) {
		out[[i]] <- mySplit(strVec, keepList[[i]])		
	}
	out <- .Call( "myPasteWithSep", as.character(unlist(out)), length(out), coll, PACKAGE = "sdcTable")
}
# mySplitIndicesList("112233444", list(1:3, 5:6, 7:8))

# check ob 'x' ganzzahlig ist
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  {
	abs(x - round(x)) < tol
}

is.zero <- function(x, tol = .Machine$double.eps^0.5) {
	abs(x - 0) < tol
}

is.one <- function(x, tol = .Machine$double.eps^0.5) {
	abs(x - 1) < tol
}

# welche Variable soll als Branching_Variable verwendet werden?
getBranchingVariable <- function(sol, alreadyBranched, primSupps) {
	ind <- setdiff(1:length(sol), c(alreadyBranched, primSupps))	
	branchVar <- ind[which.min(0.5 - sol[ind])]
	branchVar
}

my.Rglpk_solve_LP <- function(obj, mat, dir, rhs, types = NULL, max = FALSE, bounds = NULL, verbose = FALSE) {
	if (!identical(max, TRUE) && !identical(max, FALSE)) 
		stop("'Argument 'max' must be either TRUE or FALSE.")
	direction_of_optimization <- as.integer(max)
	if (!identical(verbose, TRUE) && !identical(verbose, FALSE)) 
		stop("'Argument 'verbose' must be either TRUE or FALSE.")
	if ( !class(mat) == "simpleTriplet" ) {
		stop("mat must be of class 'simpleTriplet'")
	}
	
	verb <- as.integer(verbose)
	n_of_constraints <- length(dir)
	direction_of_constraints <- match(dir, c("<", "<=", ">", ">=", "=="))
	if (any(is.na(direction_of_constraints))) 
		stop("Argument 'dir' must be either '<', '<=', '>', '>=' or '=='.")
	n_of_objective_vars <- length(obj)
	if (is.null(types)) 
		types <- "C"
	if (any(is.na(match(types, c("I", "B", "C"), nomatch = NA)))) 
		stop("'types' must be either 'B', 'C' or 'I'.")
	types <- rep(types, length.out = n_of_objective_vars)
	integers <- types == "I"
	binaries <- types == "B"
	is_integer <- any(binaries | integers)
	bounds <- as.glp_bounds(as.list(bounds), n_of_objective_vars)
	x <- glp_call_interface(
		obj, 
		n_of_objective_vars, 
		sdcTable:::get.simpleTriplet(mat, type='rowInd', input=list()), 			
		sdcTable:::get.simpleTriplet(mat, type='colInd', input=list()), 
		sdcTable:::get.simpleTriplet(mat, type='values', input=list()), 
		length(sdcTable:::get.simpleTriplet(mat, type='values', input=list())), 		
		
		rhs, direction_of_constraints, n_of_constraints, is_integer, 
		integers, binaries, 
		direction_of_optimization, 
		bounds[,1L], 
		bounds[,2L],
		bounds[,3L], 
		verb)
	solution <- x$lp_objective_vars_values
	solution[integers | binaries] <- round(solution[integers | binaries])
	status <- as.integer(x$lp_status != 5L)
	list(optimum = sum(solution * obj), solution = solution, status = status)
}
environment(my.Rglpk_solve_LP) <- environment(Rglpk_solve_LP) 

# calculates years, weeks, days, hours, minutes and seconds from integer number 
# secs: count of elapsed seconds (proc.time()[3])
# returns also a formatted string
formatTime <- function(secs){
	time.vec <- rep(NA, 6)
	names(time.vec) <- c('seconds', 'minutes','hours', 'days','weeks','years')
	
	secs <- ceiling(secs)
	
	time.vec['years'] <- floor(secs / (60*60*24*7*52))
	if ( time.vec['years'] > 0 ) {
		secs <- secs - (time.vec['years'] * (60*60*24*7*52))
	}
	
	time.vec['weeks'] <- floor(secs / (60*60*24*7))	
	if ( time.vec['weeks'] > 0 ) {
		secs <- secs - (time.vec['weeks'] * (60*60*24*7))
	}	

	time.vec['days'] <- floor(secs / (60*60*24))
	if ( time.vec['days'] > 0 ) {
		secs <- secs - time.vec['days']*(60*60*24)
	}
	
	time.vec['hours'] <- floor(secs / (60*60))
	if ( time.vec['hours'] > 0 ) {
		secs <- secs - time.vec['hours']*(60*60)	
	}
	
	time.vec['minutes'] <- floor(secs / (60))
	if ( time.vec['minutes'] > 0 ) {
		secs <- secs - time.vec['minutes']*(60)
	}
	
	time.vec['seconds'] <- secs
	time.vec <- rev(time.vec)

	# time str #
	x <- time.vec[time.vec!=0]
	shortNames <- sapply(1:length(x), function(y) { substr(names(x)[y], 1, nchar(names(x)[y])-1)  } )
	
	time.str <- NULL
	for ( i in seq_along(names(x))) {
		
		if ( length(x) == 1 ) {
			if ( x[i] > 1 ) {
				time.str <- paste(time.str, x[i], " ", names(x[i]), sep="")	
			} else {
				time.str <- paste(time.str, x[i], " ", shortNames[i], sep="")	
			}			
		} 
		else {
			
			if ( names(x)[i]=="seconds") {
				if ( x[i] > 1 ) {
					time.str <- paste(time.str, "and", x[i], names(x[i]), sep=" ")	
				} else {
					time.str <- paste(time.str, "and", x[i], shortNames[i], sep=" ")	
				}
				
			} else {
				if ( x[i] > 1 ) {
					time.str <- paste(time.str, x[i], " ", names(x[i]), sep="")	
				} else {
					time.str <- paste(time.str, x[i], " ", shortNames[i], sep="")	
				}
				
				if ( i != length(x)-1 ) {
					time.str <- paste(time.str,", ", sep="")	
				}
			}			
		}
	}
	return(list(time.vec=time.vec, time.str=time.str))
}

# create default parameter objects suitable for primary|secondary suppression
# if selection == 'control.primary': set arguments suitable for primary suppression
# if selection == 'control.secondary': set arguments suitable for secondary suppression
genParaObj <- function(selection, ...) {
	controlPrimary <- function(...) {		
		### setDefaults ###
		paraObj <- list()
		
		# freq.rule
		paraObj$maxN <- 3
		paraObj$allowZeros <- FALSE
		
		# p-percent rule
		paraObj$p <- 80		
		
		# n,k rule
		paraObj$n <- 2
		paraObj$k <- 85
		
		paraObj$numVarInd <- NA
		
		newPara <- list(...)
		
		indexNumVarIndices <- which(names(newPara) == "numVarIndices")
		if ( length(indexNumVarIndices) == 0 ) {
			stop("genPara (type=='control.primary'): parameter 'numVarIndices' must be specified\n")
		} else {
			numVarIndices <- newPara[[indexNumVarIndices]]	
		}
		
		for ( i in seq_along(newPara) ) {
			m <- match(names(newPara)[i], names(paraObj))
			if ( !is.na(m) ) {
				paraObj[[m]] <- newPara[[i]]
			}
		}
		
		if ( any(sapply(paraObj, length)!=1) ) {
			stop("genPara (type=='control.primary'): arguments for primary suppression are not valid!\n")
		}
		if ( !is.logical(paraObj$allowZeros) ) {
			stop("genPara (type=='control.primary'): argument 'allowZeros' must be logical!\n")
		}
		
		if ( !all(c(is.numeric(paraObj$maxN), is.numeric(paraObj$p), is.numeric(paraObj$n), is.numeric(paraObj$k))) ) {
			stop("genPara (type=='control.primary'): arguments 'maxN', 'p', 'n' and 'k' must be numeric!\n")
		}
		
		if ( paraObj$k < 1 | paraObj$k >= 100) {
			stop("genPara (type=='control.primary'):argument 'k' must be >= 1 and < 100!\n")
		}			
		if ( paraObj$p < 1 | paraObj$p >= 100) {
			stop("genPara (type=='control.primary'):argument p must be >= 1 and < 100!\n")
		}	
		if ( !is.na(paraObj$numVarInd) ) {
			if ( !paraObj$numVarInd %in% 1:length(numVarIndices) ) {
				stop("genPara (type=='control.primary'):argument 'numVarInd' must be >= 1 and <=",length(numVarIndices),"!\n")
			}				
		}
		return(paraObj)
	}
	
	### create a parameter list with (...) changing the default-values -> used in protectTable()
	controlSecondary <- function(...) {
		### setDefaults ###
		paraObj <- list()
		
		# general parameter
		paraObj$method <- NA
		paraObj$verbose <- FALSE
		paraObj$save <- FALSE		
		paraObj$solver <- "glpk"
		
		# HITAS|OPT - parameter
		paraObj$maxIter <- 10
		paraObj$timeLimit <- NULL 
		paraObj$maxVars <- NULL 
		paraObj$fastSolution <- FALSE 
		paraObj$fixVariables <- TRUE 
		paraObj$approxPerc <- 10
		
		# HYPERCUBE - parameter
		paraObj$protectionLevel <- 80 
		paraObj$suppMethod <- "minSupps" 
		paraObj$suppAdditionalQuader <- FALSE
		
		# protectLinkedTables
		paraObj$maxIter <- 5
		
		newPara <- list(...)
		for ( i in seq_along(newPara) ) {
			m <- match(names(newPara)[i], names(paraObj))
			if ( !is.na(m) ) {
				paraObj[[m]] <- newPara[[i]]
			}
		}
		
		### checks
		if ( any(sapply(paraObj, length)!=1) ) {
			stop("genPara (type=='control.secondary'): arguments controlObj for sdc-procedure are not valid!\n")
		}	
		if ( !all(c(is.numeric(paraObj$maxIter), is.numeric(paraObj$approxPerc), is.numeric(paraObj$protectionLevel), is.numeric(paraObj$maxIter))) ) {
			stop("genPara (type=='control.secondary'): arguments 'maxIter', 'maxIter', 'protectionLevel' and 'maxIter' must be numeric!\n")
		}
		if ( !all(c(is.logical(paraObj$verbose), is.logical(paraObj$save), is.logical(paraObj$fastSolution), is.logical(paraObj$fixVariables), is.logical(paraObj$suppAdditionalQuader))) ) {
			stop("genPara (type=='control.secondary'): arguments 'verbose', 'save', 'fastSolution' 'fixVariables' and 'suppAdditionalQuader' must be numeric!\n")
		}	
		if ( !is.null(paraObj$timeLimit) && !paraObj$timeLimit %in% 1:3000 ) {
			stop("genPara (type=='control.secondary'): argument 'timeLimit' must be >= 1 and <= 3000 minutes!\n")
		}
		if ( !length(paraObj$approxPerc) & !paraObj$approxPerc %in% 1:100 ) {
			stop("genPara (type=='control.secondary'): argument 'approxPerc' must be >= 1 and <= 100!\n")
		}
		if ( !paraObj$method %in% c('HITAS', 'HYPERCUBE', 'OPT') ) {
			stop("genPara (type=='control.secondary'): 'method' must be either 'HITAS', 'HYPERCUBE' or 'OPT'!\n")
		}
		if ( !paraObj$suppMethod %in% c('minSupps', 'minSum', 'minSumLogs') ) {
			stop("genPara (type=='control.secondary'): 'suppMethod' must be either 'minSupps', 'minSum' or 'minSumLogs'!\n")
		}
		return(paraObj)
	}
	
	if ( !selection %in% c('control.primary', 'control.secondary') ) {
		stop("genPara:: argument 'selection' must be either 'control.primary' or 'control.secondary'!\n")
	}
	
	if ( selection == 'control.primary' ) {
		paraObj <- controlPrimary(...)	
	}
	if ( selection == 'control.secondary' ) {
		paraObj <- controlSecondary(...)	
	}		
	return(paraObj)
}