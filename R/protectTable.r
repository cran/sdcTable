protectTable <- function(outObj, method, ...) {
	# if levelObj is not NULL -> complete problem
	# else problem for subTable only
	f.genSubTab <- function(fullTabObj, strIDs, levelObj=NULL) {
		subTab <- list()			
		if ( is.null(levelObj) ) {
			subTab$strID <- strIDs
			subTabInd <- match(subTab$strID, fullTabObj$strID)
		}
		else {	
			subTabInd <- 1:length(fullTabObj$strID)
			subTab$strID <- fullTabObj$strID[subTabInd]
		} 
		
		subTab$Freq <- fullTabObj$Freq[subTabInd]
		subTab$w <- fullTabObj$w[subTabInd]
		subTab$numVal <- fullTabObj$numVal[subTabInd]
		subTab$status <- fullTabObj$status[subTabInd]
		subTab$lb <- fullTabObj$lb[subTabInd]
		subTab$ub <- fullTabObj$ub[subTabInd]
		subTab$LPL <- fullTabObj$LPL[subTabInd]
		subTab$UPL <- fullTabObj$UPL[subTabInd]
		subTab$SPL <- fullTabObj$SPL[subTabInd]		
		subTab$UB <- fullTabObj$UB[subTabInd]
		subTab$LB <- fullTabObj$LB[subTabInd]
		
		# TODO: protectTable() and HITAS() needs to be expanded to take care of this
		#subTab$weight <- "default"
		return(subTab)		
	}
	
	# calculate the information neccessary to perform HITAS|HYPERCUBE
	# this functions generates the correct order of the subtables to be 
	# processed by HITAS|HYERCUBE to generate valid suppression patterns
	# Since it is used for both approaches, it was separated from the
	# main code and f.getSplitTables() was introduced
	f.genSplitTables <- function(fullTabObj, levelObj) {
		fn.tabLevels <- function(levs, levelObj) {
			fn.prepareResult <- function(result) {
				x <- list()
				for ( i in 1:length(result))
					x[[i]] <- 1:length(result[[i]])
				combs <- expand.grid(x)
				
				vecs <- list()
				for ( i in 1:nrow(combs)) {
					str <- paste("vecs[[",i,"]] <- expand.grid(result[[1]][[",combs[i,1],"]]", sep="")
					if ( ncol(combs) > 1) {
						for ( j in 2:ncol(combs) ) {
							str <- paste(str, ", result[[",j,"]][[",combs[i,j],"]]", sep="")
						}	
					}
					str <- paste(str, ", sep='')", sep="")		
					eval(parse(text=str))
					vecs[[i]] <- apply(vecs[[i]], 1, paste, collapse="")
					
				}
				vecs
			}			
			result <- list()
			for ( i in 1:length(levelObj) ) {
				result[[i]] <- list()
				index <- 1:length(levelObj[[i]]$levelsOrig)
				index <- index[which(levelObj[[i]]$levelsOrig %in% c(levs[i], levs[i]-1))]
				levOrig <- levelObj[[i]]$levelsOrig[index]
				codesStandard <- levelObj[[i]]$codesStandard[index]
				
				
				if ( levs[i] == 1) {
					result[[i]] <- codesStandard
				}
				else {
					diffs <- c(0,diff(levOrig))
					checkInd <- which(diffs == 1)-1
					out <- data.frame(index=index, levOrig=levOrig, codesStandard=codesStandard, ind=NA)
					out$ind[checkInd] <- 1
					
					checkInd <- c(checkInd, length(index))
					splitVec <- rep(0, length(index))
					for ( z in 2:length(checkInd) ) {
						if ( z < length(checkInd) )
							splitVec[checkInd[z-1]:(checkInd[z]-1)] <- z-1
						else
							splitVec[checkInd[z-1]:(checkInd[z])] <- z-1
					}
					spl <- split(index, splitVec)
					
					counter <- 1
					for ( z in 1:length(spl) ) {
						rowInd <- match(spl[[z]], out$index)
						tmp <- out[rowInd,]
						if ( any(tmp[,"levOrig"]==levs[i]) ) {					
							tmp <- tmp[1:(max(which(tmp$levOrig==levs[i]))),]
							result[[i]][[counter]] <- sort(unique(as.character(tmp$codesStandard)))
							counter <- counter + 1	
						}
					}				
				}
			}
			result	<- fn.prepareResult(result)
			result
		}
		
		# 1) create classes and groups
		tmpDat <- expand.grid(lapply(levelObj, function(x) { 1:length(x$levelStructure) } ))
		groups <- apply(tmpDat, 1, function(x) { paste(x, collapse="-")})
		classes <- apply(tmpDat, 1, sum)
		sortOrder <- order(classes)
		classes <- classes[sortOrder]
		classesUnique <- unique(classes)
		groups <- groups[sortOrder]
		splitGroups <- split(groups, classes)
		
		# 2) create tables for all classes and groups (fn.tabLevels)
		counter <- 0
		out <- list()
		out$groups <- groups
		out$indices <- list()
		for ( i in 1:length(groups) ) {
			out$indices[[i]] <- list()
			levs <- as.integer(unlist(sapply(groups[i], strsplit, "-")))		
			xx <- fn.tabLevels(levs, levelObj)	
			for ( j in 1:length(xx) ) {
				out$indices[[i]][[j]] <- xx[[j]]	
				counter <- counter + 1
			}
			
		}
		out$nrTables <- counter	
		out$nrGroups <- length(groups)
		out
	}	
	
	# finalize output objects
	# -> use objects from protectTable()
	# -> pretty formatting [TODO]
	# -> options on how to export the data [TODO]
	f.finalizeOutput <- function(fullTabObj, strObj, levelObj, method, solver, addSuppsInfo, end, start) {	
		strInfo <- strObj$strInfo
		
		nrPrimarySupps <- length(which(fullTabObj$status == "u"))
		nrSecondSupps <- length(which(fullTabObj$status == "x"))
		
		allCodes <- list()
		for ( i in 1:length(levelObj) ) {
			allCodes[[i]] <- list()
			if ( length(levelObj[[i]]$dups) == 0 ) {
				x <- levelObj[[i]]$codesStandard[match(levelObj[[i]]$dupsUp, levelObj[[i]]$codesOrig)]
				levDefault <- c(levelObj[[i]]$codesStandard, as.character(sapply(x, calcUpperLev, levelObj[[i]])))
				codesOrig <- c(levelObj[[i]]$codesOrig, levelObj[[i]]$dups)	
				allCodes[[i]]$codesStandard <- levDefault
				allCodes[[i]]$codesOrig <- codesOrig
			}
			else {
				dupsUp <- levelObj[[i]]$dupsUp
				dups <- levelObj[[i]]$dups
				codesOrig <- levelObj[[i]]$codesOrig
				levDefault <- levelObj[[i]]$codesStandard
				while( length(dupsUp) > 0 ) {
					matchInd <- match(dupsUp, codesOrig)
					x <- levDefault[na.omit(matchInd)]
					
					levDefault <- c(levDefault, as.character(sapply(x, calcUpperLev, levelObj[[i]])))
					codesOrig <- c(codesOrig, dups[!is.na(matchInd)])	

					dupsUp <- setdiff(dupsUp,  dupsUp[!is.na(matchInd)])
					dups <- setdiff(dups,  dups[!is.na(matchInd)])
				}	
				allCodes[[i]]$codesStandard <- levDefault
				allCodes[[i]]$codesOrig <- codesOrig
			}			
		}	
		
		outObj <- fullTabObj
		for ( i in 1:length(levelObj) ) {
			if ( length(levelObj[[i]]$dups) > 0 ) {
				orderInd <- order(match(levelObj[[i]]$dupsUp, levelObj[[i]]$codesStandard), decreasing=TRUE)
				levelObj[[i]]$dups <- levelObj[[i]]$dups[orderInd]
				levelObj[[i]]$dupsUp <- levelObj[[i]]$dupsUp[orderInd]

				for ( j in 1:length(levelObj[[i]]$dups) ) {
					code <- levelObj[[i]]$dupsUp[j]
					codeS <- allCodes[[i]]$codesStandard[match(code, allCodes[[i]]$codesOrig)]
					codeU <- calcUpperLev(codeS, levelObj[[i]])
					
					addIndex <- which(substr(outObj$strID, strInfo[[i]][1], strInfo[[i]][2]) == codeS)
									
					if ( length(addIndex) == 0 ) 
						stop("Error: i=",i,"; j=",j,"!\n")

					strIDNew <- outObj$strID[addIndex]
					substr(strIDNew, strObj$strInfo[[i]][1], strObj$strInfo[[i]][2]) <- rep(codeU, length(addIndex))
					outObj$strID <- c(outObj$strID, strIDNew)
					outObj$Freq <- c(outObj$Freq, outObj$Freq[addIndex])
					outObj$w <- c(outObj$w, outObj$Freq[addIndex])		
					outObj$status <- c(outObj$status, outObj$status[addIndex])	
					outObj$numVal <- c(outObj$numVal, outObj$numVal[addIndex])
					outObj$lb <- c(outObj$lb, outObj$lb[addIndex])
					outObj$ub <- c(outObj$ub, outObj$ub[addIndex])
					outObj$LPL <- c(outObj$LPL, outObj$LPL[addIndex])
					outObj$UPL <- c(outObj$UPL, outObj$UPL[addIndex])
					outObj$SPL <- c(outObj$SPL, outObj$SPL[addIndex])
					outObj$UB <- c(outObj$UB, outObj$UB[addIndex])
					outObj$LB <- c(outObj$LB, outObj$LB[addIndex])
				}
			}			
		}
		
		# Merge Bezeichnungen und Codes
		outCodes <- list()
		for ( i in 1:length(levelObj) ) {
			outCodes[[i]] <- substr(outObj$strID, strInfo[[i]][1], strInfo[[i]][2])
			outCodes[[i]] <- allCodes[[i]]$codesOrig[match(outCodes[[i]], allCodes[[i]]$codesStandard)]
		}
		names(outCodes) <- unlist(lapply(levelObj, function(x) x$varName))
		outObj <- c(outObj, outCodes)
		
		if ( all(outObj$numVal %in% c(0, NA)) )
			outObj$numVal <- rep(NA, length(outObj$numVal))
		
		protectedObj <- list()
		protectedObj$time <- formatC(round((end-start)/60,2), digits=2, format="f")
		protectedObj$method <- method
		protectedObj$nrPrimarySupps <- nrPrimarySupps
		protectedObj$nrSecondSupps <- nrSecondSupps
		protectedObj$outObj <- outObj
		protectedObj$suppsInfo <- addSuppsInfo
		protectedObj$solver <- solver
		protectedObj$levelObj <- levelObj
		protectedObj$strObj <- strObj
		class(protectedObj) <- "safeTable"
		return(protectedObj)				
		protectedObj	
	}
	
	#########################
	### Start the program ###
	#########################
	
	### start: preparations ###	
	# adding additional specific parameters for HITAS|GHMITER
	suppMethod = list(...)$suppMethod
	protectionLevel = list(...)$protectionLevel
	allowZeros = list(...)$allowZeros
	randomResult = list(...)$randomResult
	solver = list(...)$solver
	cpus = list(...)$cpus
	if( !method %in% c("HYPERCUBE", "HITAS", "OPT") ) 
		stop("===> Error: please choose a valid method for secondary cell suppression! Choices are \"HYPERCUBE\" and \"HITAS\"!\n"); flush.console()
	
	start <- as.numeric(proc.time())[3]
	strObj <- outObj$strObj
	strInfo <- strObj$strInfo
	levelObj <- outObj$levelObj	
	fullTabObj <- outObj$fullTabObj			
	
	if ( method %in% c("HITAS", "OPT") ) {
		if ( !is.null(solver) ) {
			if ( !solver %in% c("lpsolve", "glpk", "symphony" ,"cplex") ) {
				cat("===> Warning: the solver you have chosen must be one of \"lpsolve\"|\"symphony\"|\"cplex\"|\"glpk\"! Therefore, the default choice \"glpk\" will be used.\n"); flush.console()					
				solver <- "glpk"
			}
		}
		if( is.null(solver) )  {
			cat("===> Note: the default choice \"glpk\" for parameter \"solver\" will be used!\n"); flush.console()
			solver <- "glpk"
		}				
		if ( solver == "cplex" &  get("indCplex", pos=which(search()=="myGlobalEnv")) == FALSE ) {
			solver <- "gplk"
			cat("===> Note: the default choice \"glpk\" for parameter \"solver\" will be used because library 'Rcplex' is not available!\n"); flush.console()
		}
		if ( solver == "symphony" & get("indSymphony", pos=which(search()=="myGlobalEnv")) == FALSE ) {
			solver <- "gplk"
			cat("===> Note: the default choice \"glpk\" for parameter \"solver\" will be used because library 'RSymphony' is not available!\n"); flush.console()
		}			
		if ( solver == "lpsolve" & get("indLpSolveAPI", pos=which(search()=="myGlobalEnv")) == FALSE ) {
			solver <- "gplk"
			cat("===> Note: the default choice \"glpk\" for parameter \"solver\" will be used because library 'lpSolveAPI' is not available!\n"); flush.console()
		}				
		# check if fullTabObj is a valid object
		if ( any(
				is.null(fullTabObj$lb), 
				is.null(fullTabObj$ub), 
				is.null(fullTabObj$LPL), 
				is.null(fullTabObj$UPL)) )
			stop("Error: please add information about lower|upper bounds|protection levels to fullTabObj using setBounds()"); flush.console()
	}	
	if( is.null(cpus) )  {
		cat("===> Note: the default choice \"1\" for parameter \"cpus\" will be used!\n"); flush.console()
		cpus <- 1
	}			
	if ( cpus > 1 & get("indSnowfall", pos=which(search()=="myGlobalEnv")) == FALSE ) {
		cpus <- 1	
		cat("===> Note: the default choice \"1\" for parameter \"cpus\" will be used because library 'snowfall' is not availabe!\n"); flush.console()
	}		
	
	if ( method == "HYPERCUBE" ) {
		solver <- NULL
		if( is.null(suppMethod) )
			suppMethod <- "minSupps"; cat("===> NOTE: the default value \"minSupps\" for parameter \"suppMethod\" will be used!\n")
		if( is.null(protectionLevel) )
			protectionLevel <- 80; cat("===> NOTE: the default value of 80 for parameter \"protectionLevel\" will be used!\n")
		if( is.null(allowZeros) )	
			allowZeros <- FALSE; cat("===> NOTE: the default value \"FALSE\" for parameter \"allowZeros\" will be used!\n")
		if( is.null(randomResult) )	
			randomResult <- FALSE; cat("===> NOTE: the default value \"FALSE\" for parameter \"randomResult\" will be used!\n")
	}			
	### end: preparations ###
	
	### start: declare common objects ###
	# save original primary suppressed cells
	origPrimSupps <- which(fullTabObj$status=="u")
	origSecondSupps <- which(fullTabObj$status=="x")
	origProtectedCells <- which(fullTabObj$status=="z")			
	addSuppsInfo <- suppPattern <- NA		
	### end: declare common objects ###
	
	if ( any(fullTabObj$status == "u") ) {
		addSuppsInfo <- suppPattern <- list()
		
		### start: HITAS procedure ###
		if ( method == "HITAS" ) {
			suppPattern <- addSuppsInfo <- list()
			addSupps <- c()
			splitInfo <- f.genSplitTables(fullTabObj, levelObj)
			nrGroups <- splitInfo$nrGroups
			counter <- 0
			i <- j <- 1		
			while ( i <= nrGroups ) {
				cat("working on group",i,"|",nrGroups,"\n")
				j <- 1
				lenTabs <- length(splitInfo$indices[[i]])
				while ( j <= lenTabs ) {	
					error <- FALSE
					counter <- counter + 1 
					addSuppsInfo[[counter]] <- list()
					subTab <- f.genSubTab(fullTabObj, splitInfo$indices[[i]][[j]], levelObj=NULL)
					# do primary supps (status=="u") exist in marginals? if so -> "x"
					cellInfo <- isMarginalSum(subTab$strID, strInfo) 
					indX <- which(subTab$status=="x" & subTab$Freq != 1)
					indXTot <- indX[which(indX %in% cellInfo$indexTotCells)]					
					if ( length(indXTot) > 0 )
						subTab$status[indXTot] <- "u"
					
					# standard case: primary suppressions exist
					if ( any(subTab$status %in% c("u","x")) && any(subTab$status %in% c("z","s")) ) {						
						M <- genMatM(subTab$strID, strObj$strInfo)
						subTab <- protectHITAS(subTab, M, solver, verbose=FALSE, cpus, method)
						
						#cat("temporarily saving protected subTab....")
						#save(subTab, file=paste("protectedSubTab-",i,"-",j,".RData", sep=""))
						#cat("Done!\n")						
						
						if ( is.null(subTab) ) {
							# subTab was NULL -> recalculate 
							subTab <- f.genSubTab(fullTabObj, splitInfo$indices[[i]][[j]], levelObj=NULL)
							totCellIndexChanged <- cellInfo$indexTotCells[na.omit(match(which(subTab$status=="z"), cellInfo$indexTotCells))]
							subTab$status[totCellIndexChanged] <- "s"
							
							# protect the "relaxed" subTab
							subTab <- protectHITAS(subTab, M, solver, verbose=FALSE, cpus, method)
							
							# if there is a solution now
							if ( !is.null(subTab) ) {
								changed <- subTab$strID[totCellIndexChanged[subTab$status[totCellIndexChanged]=="x"]]
								if ( length(changed) > 0 ) {
									addSupps <- c(addSupps, match(changed, fullTabObj$strID))
									fullTabObj <- outObj$fullTabObj	
									fullTabObj$status[addSupps] <- "u"
									error <- TRUE
								}									
							}	
							else
								stop("something is terribly wrong!\n")
						}
					}				
					if ( error == TRUE ) {
						j <- lenTabs+1
						i <- 0
					}
					else {
						j <- j + 1						
						# else : we update the suppression pattern
						secondSuppInd <- match(subTab$strID[which(subTab$status=="x")], fullTabObj$strID)
						addSuppsInfo[[counter]]$primSupps <- NA
						addSuppsInfo[[counter]]$secondSupps <- NA	
						if ( length(secondSuppInd) > 0 ) {
							fullTabObj$status[secondSuppInd] <- "x"
							addSuppsInfo[[counter]]$primSupps <- subTab$strID[subTab$status=="u"]
							addSuppsInfo[[counter]]$secondSupps <- fullTabObj$strID[secondSuppInd]					
						}
						
						# new in HITAS: 
						# all levels already dealt with need to be set to 'publishable' temporarily
						vecNotChangable <- match(subTab$strID[which(subTab$status%in%c("s","z"))], fullTabObj$strID)
						if ( length(vecNotChangable) > 0 )
							fullTabObj$status[vecNotChangable] <- "z"	
					}				
				}		
				i <- i + 1
			}			
		}	
		### end: HITAS procedure ###
		
		
		### start: HYPERCUBE procedure ###
		if ( method == "HYPERCUBE" ) {
			runInd <- TRUE
			nrRuns <- 1
			while ( runInd == TRUE ) { 	
				primSupps <- which(fullTabObj$status=="u")
				secondSupps <- which(fullTabObj$status=="x")
				protectedCells <- which(fullTabObj$status=="z")
				if ( any(fullTabObj$status == "u") ) {
					suppPattern <- addSuppsInfo <- list()
					splitInfo <- f.genSplitTables(fullTabObj, levelObj)
					nrGroups <- splitInfo$nrGroups
					counter <- 0
					i <- j <- 1		
					for ( i in 1:nrGroups ) {
						lenTabs <- length(splitInfo$indices[[i]])
						for ( j in 1:lenTabs ) {	
							error <- FALSE
							counter <- counter + 1 
							addSuppsInfo[[counter]] <- list()
							addSuppsInfo[[counter]]$primSupps <- NA
							addSuppsInfo[[counter]]$secondSupps <- NA	
							
							subTab <- f.genSubTab(fullTabObj, splitInfo$indices[[i]][[j]], levelObj=NULL)
							
							# do primary supps (status=="u") exist in marginals? if so -> "x"
							cellInfo <- isMarginalSum(subTab$strID, strInfo) 
							
							indX <- which(subTab$status=="x" & subTab$Freq != 1)
							indXTot <- indX[which(indX %in% cellInfo$indexTotCells)]
							
							if ( length(indXTot) > 0 ) 
								subTab$status[indXTot] <- "u"							
							# standard case: primary suppressions exist
	
							if ( any(subTab$status =="u") && any(subTab$status != "u") ) {
								subTab <- protectHYPERCUBE(subTab, levelObj, strObj, allowZeros, protectionLevel, suppMethod, debug=FALSE)
								if ( !is.null(subTab) )  {
									secondSuppInd <- match(subTab$strID[which(subTab$status=="x")], fullTabObj$strID)
									if ( length(secondSuppInd) > 0 ) {
										fullTabObj$status[secondSuppInd] <- "x"
										addSuppsInfo[[counter]]$primSupps <- subTab$strID[subTab$status=="u"]
										addSuppsInfo[[counter]]$secondSupps <- fullTabObj$strID[secondSuppInd]					
									}
								} else 
									stop("Error with subTab (i=",i,"|j=",j,")!\n")
							}			
						}		
					}		
				}
				# else: return input-data as they are, nothing to do
				end <- as.numeric(proc.time())[3]	
				
				newSupps <- setdiff(which(fullTabObj$status %in% c("u","x")),  c(primSupps, secondSupps))
				
				if ( length(newSupps) == 0 ) {
					runInd <- FALSE
				}
				else {
					protectionLevel <- 1
					fullTabObj$status <- rep("s", length(fullTabObj$status))
					fullTabObj$status[newSupps] <- "u"
					fullTabObj$status[unique(c(primSupps, secondSupps))] <- "x"
					fullTabObj$status[protectedCells] <- "z"
					nrRuns <- nrRuns + 1
				}
			}			
			
			# restore original primary suppressed cells
			# index of cells that were set to "u" in the heuristic process
			# but are in fact secondary suppressions
			index <- which(fullTabObj$status %in% c("x","u"))
			
			# restore original state of cells that are forced to be published
			fullTabObj$status <- rep("s", length(fullTabObj$status))
			if ( length(origProtectedCells) > 0 )
				fullTabObj$status[origProtectedCells] <- "z"	
			
			if ( length(index) > 0 )
				fullTabObj$status[index] <- "x"
			
			if ( length(origPrimSupps) > 0 )
				fullTabObj$status[origPrimSupps] <- "u"			
			
		}		
		### end: HYPERCUBE procedure ###	
		
		
		### start: OPT procedure ###
		if ( method == "OPT" ) {
			#subTab <- f.genSubtab(fullTabObj, currentLevel=NULL, dimObj=NULL, levelObj)
			subTab <- f.genSubTab(fullTabObj, strIDs=NULL, levelObj)
			M <- genMatMFull(subTab$strID, levelObj)			
			subTab <- protectHITAS(subTab, M, solver, verbose=FALSE, cpus, method)
			secondSuppInd <- match(subTab$strID[which(subTab$status=="x")], fullTabObj$strID)
			if ( length(secondSuppInd) > 0 ) {
				fullTabObj$status[secondSuppInd] <- "x"
				addSuppsInfo$primSupps <- subTab$strID[subTab$status=="u"]
				addSuppsInfo$secondSupps <- fullTabObj$strID[secondSuppInd]					
			}
			else {
				addSuppsInfo$primSupps <- NA
				addSuppsInfo$secondSupps <- NA							
			}			
		}			
		### end: OPT procedure ###				
		
		### start: restore original cell states ###
		# restore original primary suppressed cells
		# index of cells that were set to "u" in the heuristic process
		# but are in fact secondary suppressions
		index <- which(fullTabObj$status %in% c("x","u"))
		
		# restore original state of cells that are forced to be published
		fullTabObj$status <- rep("s", length(fullTabObj$status))
		if ( length(origProtectedCells) > 0 )
			fullTabObj$status[origProtectedCells] <- "z"	
		
		if ( length(index) > 0 )
			fullTabObj$status[index] <- "x"
		
		if ( length(origPrimSupps) > 0 )
			fullTabObj$status[origPrimSupps] <- "u"		
		### end: restore original cell states ###
		
	}
	# else: return input-data as they are, nothing to do
	end <- as.numeric(proc.time())[3]
	
	# TODO: what should we do with NA-Cells? 
	# possibility: mark them as secondary suppressions with temporary value > 1	
	protectedObj <- f.finalizeOutput(fullTabObj, strObj, levelObj, method, solver, addSuppsInfo, end, start)
	cat("\n===>",protectedObj$nrPrimarySupps,"primary sensitive cells have been protected with",protectedObj$nrSecondSupps,"secondary cell-suppressions using", method, "algorithm. [Finished]\n")			
	return(protectedObj)	
}
