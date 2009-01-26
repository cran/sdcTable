processTableHYPERCUBE <- function(fullData, allowZeros=TRUE, randomResult=FALSE, suppMethod="minSupps", protectionLevel=80) {
	recodeIndexVars <- function(subtab, indexvars) {
		tmp <- apply(subtab[,indexvars], 2, function(x) { sort(unique(x)) } )
		nc <- ncol(subtab)
		subtab2 <- cbind(subtab, matrix(0, nrow=nrow(subtab), ncol=length(indexvars)))
		
		if(is.matrix(tmp)) {
			for (i in 1:length(indexvars)) {
				for (z in 1:nrow(tmp)) {
					subtab2[which(as.character(subtab[,indexvars[i]]) == tmp[z,i]), nc+i] <- z
				}
			}
		}
		else {
			for (i in 1:length(indexvars)) {
				for (z in 1:length(tmp[[i]])) {
					subtab2[which(subtab[,indexvars[i]] == tmp[[i]][[z]]), nc+i] <- z
				}
			}
		}
		subtab2
	}	
	
	# safeQuaderSubtable protects all primary/secondary suppressed values within a subtable as well as single cells according to the specified suppression method and protection level (protectionLevel)
	## subtab: the current subtable from fullData$data with suppressed values to protect
	## allowZeros: TRUE/FALSE - are empty cells allowed in the suppression scheme
	## randomResult: TRUE/FALSE - if several suppression schemes are equal, should one pick a random scheme or always the first one?
	## suppMethod: "minSupps"|"minVals" - choose suppression scheme with minimum (additional) suppressions or with minimum (additional) amount of information
	## protectionLevel: protection level 
	algorithmGHMITER <- function(subtab, allowZeros, randomResult, suppMethod, protectionLevel) {
		# for a given supressed value, calculate a list of all possible diametral indices
		diametralIndex <- function(subtab, indGeh) {
			nc <- ncol(subtab$data)
			beg <- nc - subtab$numberindexvars
			tmp <- subset(subtab$data, subtab$data[,beg+1] != indGeh[1])
			for (i in 2:subtab$numberindexvars) {
				tmp <- subset(tmp, tmp[,beg+i] != indGeh[i])
			}

			# position of the suppressed value indGeh within subtab$data
			indDat <- apply(subtab$data[,(subtab$numberindexvars+3):(ncol(subtab$data))], 1, function(x) { paste(x, collapse="") } )
			posIndGeh <- which(indDat %in% paste(indGeh, collapse=""))

			# is the value primary (P) or secondary (S) suppressed?
			SorP <- subtab$data[posIndGeh, "geh"]

			supp <- list()
			supp$indGeh <- indGeh
			supp$posIndGeh <- posIndGeh
			supp$SorP <- SorP
			supp$diametralIndices <- tmp[,(beg+1):nc]
			return(supp)
		}

		## berechnet das alte "iqs" (also die Quaderinformation für alle Quader indiziert durch g und einen diametralen Wert
		# subtab ... die Untertabelle zusätzlich Metainformation (indexvars,...)
		# supp ... output von diametralIndex(): enthält den Index des geheimen Wertes und alle diametralen Indices
	
		# calculateInformationForSuppValg() calculates all the information needed for a quader which is given by g and its diametral value
		## subtab: a subtable of fullData$data containing suppressed values
		## supp: output object of diametralIndex(). containes the index of the actual value which needs to be protected and all possible diametral indices.
		## allowZeros, protectionLevel as above
		calculateInformationForSuppValg <- function(subtab, supp, allowZeros, protectionLevel) {
		# the following functions are mostly implemented in C to speed up computing and are just auxiliary functions.
			calcQuader <- function(g, dia, numberIndexVars) {
				out <- rep(0, 2^numberIndexVars*numberIndexVars)
				q <- .C("calcQuader", g=as.integer(g),
									  dia=as.integer(dia),
									  numberIndexVars=as.integer(numberIndexVars),
									  out=as.integer(out), PACKAGE="sdcTable", NUOK=TRUE)$out
				Quader <- data.frame(matrix(q, ncol=numberIndexVars))
				return(list(Quader=Quader, QuaderVec=q))
			}

			calcNormQuader <- function(indexQuader, numberIndexVars) {
				normQuader <- .C("normQuader", indexQuader=as.integer(indexQuader),
											   numberIndexVars=as.integer(numberIndexVars),
											   lengthVec=as.integer(length(indexQuader)), PACKAGE="sdcTable", NUOK=TRUE)$indexQuader
				normQuader <- data.frame(matrix(normQuader, ncol=numberIndexVars, byrow=TRUE))
				normQuader
			}

			calcAggregationsstufen <- function(indexQuaderVec, numberIndexVars) {
				minDims <- as.integer(.C("calcMinimum", as.integer(indexQuaderVec),
														erg=as.integer(rep(0,numberIndexVars)),
														as.integer(numberIndexVars), PACKAGE="sdcTable", NUOK=TRUE)$erg)
				aggr <- .C("calcAggregationsstufen", indexQuaderVec=as.integer(indexQuaderVec),
													 minDims=as.integer(minDims),
													 numberIndexVars=as.integer(numberIndexVars),PACKAGE="sdcTable", NUOK=TRUE)$indexQuaderVec
				aggr <- matrix(aggr, ncol=numberIndexVars)
				aggr
			}

			calcIndizierung <- function(aggr, normQuader, numberIndexVars) {
				sq <- apply(aggr, 1, sum) + as.numeric(rownames(normQuader))
				indiziert <- rep("u", 2^numberIndexVars)
				indiziert[which(sq%%2==0)] <- "g"
				indiziert
			}

		   # extract indices from subtable
			extractIndicesSubtable <- function(indexQuader, lengthSub, numberIndexVars) {
				erg <- .C("extractIndicesSubtable",
									indexQuader=as.integer(indexQuader),
									lengthSub=as.integer(lengthSub),
									erg=as.integer(rep(0, numberIndexVars)),
									numberIndexVars=as.integer(numberIndexVars),
									powers=as.integer(rep(0, numberIndexVars)),
									final=as.integer(rep(0, lengthSub)), PACKAGE="sdcTable", NUOK=TRUE)
				return(list(powers=erg$powers, final=erg$final))
			}
			# extract indices from current quader
			extractIndicesAktQuader <- function(indexQuader, lengthSub, numberIndexVars, powers) {
				val <- .C("extractIndicesAktQuader",
					indexQuader=as.integer(indexQuader),
					lengthSub=as.integer(lengthSub),
					numberIndexVars=as.integer(numberIndexVars),
					powers=as.integer(powers),
					final=as.integer(rep(0, lengthSub)), PACKAGE="sdcTable", NUOK=TRUE)$final
				val
			}

			calcQuaderPosition <- function(vals, valsQ, numberIndexVars) {
				quaderPosition <- .C("calcQuaderPosition",
					as.integer(vals), as.integer(length(vals)), as.integer(valsQ),
					erg=as.integer(rep(0, 2^numberIndexVars)),
					as.integer(numberIndexVars), PACKAGE="sdcTable", NUOK=TRUE)$erg
				quaderPosition
			}

			# we create list which will contain the results for each diametral index of supp$indGeh
			ergebnis <- list()

			erg <- extractIndicesSubtable(as.vector(as.matrix(subtab$data[,(ncol(subtab$data) - subtab$numberindexvars+1):ncol(subtab$data)])), nrow(subtab$data), subtab$numberindexvars)
			vals <- erg$final
			powers <- erg$powers

			for (z in 1:nrow(supp$diametralIndices)) {
				
				# 1) we identify the quader given by supp$indGeh and supp$diametralIndices[z,]
				indexQ <- calcQuader(supp$indGeh, supp$diametralIndices[z,], subtab$numberindexvars)
				indexQuader <- indexQ$Quader
				indexQuaderVec <- indexQ$QuaderVec

				# 2) we calculate the normalized quader 
				normQuader <- calcNormQuader(indexQuaderVec, subtab$numberindexvars)

				# 3) we calculate aggregation-levels
				aggr <- calcAggregationsstufen(indexQuaderVec, subtab$numberindexvars)

				# 4) we calculate the necessary indices
				indiziert <- calcIndizierung(aggr, normQuader, subtab$numberindexvars)

				# 5) we calculate the extended indices of the current quader
				valsQ <- extractIndicesAktQuader(indexQuaderVec, 2^subtab$numberindexvars, subtab$numberindexvar, powers)

				# we calculate the position (row-indices) of the current quader within subtab$data
				quaderPosition <- calcQuaderPosition(vals, valsQ, subtab$numberindexvars)

				# 6) calculate various information about the current quader (infoQuader)
				# extract the quader from subtab$data
				aktQuader <- subtab$data[quaderPosition,]

				# which are the indices of the non-suppressed values
				indNonSupp <- which(aktQuader$geh == "")

				# 6.1) how many values need to be suppressed for this quader
				anzAddSupps <- length(indNonSupp)

				# 6.2) whats the amount of information which needs to be suppressed?
				sumAddSupps <- sum(aktQuader$val[indNonSupp])

				# 6.3) does the quader contains other single cells except for the suppressed value (supp$gehInd) to check?
				posIndGehaktQuader <- which(rownames(aktQuader) %in% rownames(subtab$data[supp$posIndGeh,]))

				# other single cells
				singleItems <- which(aktQuader$val[-posIndGehaktQuader]==1)
				if(length(singleItems) >= 1) {
					indikatorSingleItems <- TRUE
					singleItems <- aktQuader[singleItems, (subtab$numberindexvars+3):ncol(aktQuader)]
				}
				else {
					indikatorSingleItems <- FALSE
					singleItems <- NULL
				}

				# 6.4) does the quader contain empty cells?
				zeroItems <- which(aktQuader$val[-posIndGehaktQuader]==0)
				if(length(zeroItems) > 0) {
					indikatorZeroItems <- TRUE
					zeroItems <- aktQuader[zeroItems, (subtab$numberindexvars+3):ncol(aktQuader)]
				}
				else {
					indikatorZeroItems <- FALSE
					zeroItems <- NULL
				}

				# 6.5) is the quader protected enough? (protectionLevel)
				if(supp$SorP == "P") {
					if(length(which(indiziert=="u"))==0 | length(which(indiziert=="g"))==0) {
						schutz <- protectionLevel
						schutzInd <- TRUE
					}
					else {
						range <- min(aktQuader[which(indiziert=="u"),"val"]) + min(aktQuader[which(indiziert=="g"),"val"])
						X <- aktQuader[posIndGehaktQuader,"val"]
						if(X == 0) {
							if(range > min(aktQuader[which(aktQuader[,"geh"] != "P" & aktQuader[,"val"] != 0),"val"])) {
								schutz <- protectionLevel
								schutzInd <- TRUE
							}
							# !!! testing !!!
							else {
								schutz <- 0
								schutzInd <- FALSE
							}
						}
						else {
							schutz <- (100*range) / X
							ifelse(schutz >= protectionLevel, schutzInd <- TRUE, schutzInd <- FALSE)
						}
					}
				}
				# no interval protection needed for secondary suppressed values
				else {
					schutzInd <- TRUE
					schutz <- protectionLevel
				}

				# 7) return results
				if(anzAddSupps == 0 & schutzInd == TRUE & indikatorSingleItems == FALSE & indikatorZeroItems == FALSE) {
					#cat("the value is already proteced. Therefore we stop here!\n")
					return(erg = NULL)
					break
				}

				ergebnis[[z]] <- list(
										indexQuader = indexQuader,
										normQuader = normQuader,
										aggr = aggr,
										indiziert = indiziert,
										quaderPosition = quaderPosition,
										anzAddSupps=anzAddSupps,
										sumAddSupps = sumAddSupps,
										indikatorSingleItems = indikatorSingleItems,
										singleItems = singleItems,
										indikatorZeroItems = indikatorZeroItems,
										zeroItems = zeroItems,
										schutz = schutz,
										schutzInd = schutzInd
									)
			}
			return(list(subtab=subtab, supp=supp, iqsInfo=ergebnis))
		}

		# find the optimal suppression scheme
		findOptimalQuader <- function(infoQ, allowZeros, randomResult, suppMethod, protectionLevel) {
			subtab <- infoQ$subtab
			supp <- infoQ$supp
			iqsInfo <- infoQ$iqsInfo

			# Which elements of iqsInfo are NULL?
			NullElements <- which(unlist(lapply(lapply(iqsInfo, '[[', 'quaderPosition'), function(x) { length(x) } )) == 0)

			# recalc iqsInfo
			if(length(NullElements) > 0) {
				tmp <- (1:length(iqsInfo))[-NullElements]
				iqsInfo2 <- list()
				t <- 1
				for (i in tmp) {
					iqsInfo2[[t]] <- iqsInfo[[i]]
					t <- t + 1
				}
			}
			else {
				iqsInfo2 <- iqsInfo
			}

			# diametral values using alternate indices
			dia <- unlist(lapply(lapply(iqsInfo2, '[[', 'indexQuader'), function(x) { paste(x[2^length(supp$indGeh),], collapse="")} ))

			# put iqs together so that we can choose the optimal suppression scheme
			anzAddSupps <- as.numeric(as.character(do.call(rbind, lapply(iqsInfo2, '[', 'anzAddSupps'))))
			sumAddSupps <- as.numeric(as.character(do.call(rbind, lapply(iqsInfo2, '[', 'sumAddSupps'))))
			indikatorSingleItems <- as.logical(do.call(rbind, lapply(iqsInfo2, '[', 'indikatorSingleItems')))
			indikatorZeroItems <- as.logical(do.call(rbind, lapply(iqsInfo2, '[', 'indikatorZeroItems')))
			schutz <- as.numeric(as.character(do.call(rbind, lapply(iqsInfo2, '[', 'schutz'))))
			schutzInd <- as.logical(do.call(rbind, lapply(iqsInfo2, '[', 'schutzInd')))
			listInd <- 1:length(iqsInfo2)    # necessary to know which element in supps$diametralIndices defines the optimal suppression scheme

			iqs <- data.frame(listInd, dia, anzAddSupps,sumAddSupps,indikatorSingleItems,indikatorZeroItems,schutz,schutzInd)

			# all rows are removed where indicatorZeroItems != FALSE if allowZeros is FALSE
			if(allowZeros == FALSE) {
				iqs <- iqs[iqs$indikatorZeroItems == FALSE,]
			}

			# are there any suppression schemes satisfying the necessary interval protection?
			ifelse(length(which(iqs$schutzInd == TRUE)) >= 1, indexIntervallOk <- TRUE, indexIntervallOk <- FALSE)

			# do suppression schemes exist that do not contain single values? these are preferred.
			ifelse(length(which(iqs$indikatorSingleItems == FALSE)) > 0, existNonSingles <- TRUE, existNonSingles <- FALSE)
			if(existNonSingles == TRUE) {
				iqs <- iqs[iqs$indikatorSingleItems==FALSE,]
			}

			if(indexIntervallOk == TRUE) {
				iqs <- iqs[iqs$schutzInd == TRUE,]
				if(suppMethod=="minSupps") {
					iqs <- iqs[which(iqs$anzAddSupps == min(iqs$anzAddSupps)),]
					iqs <- iqs[which(iqs$sumAddSupps == min(iqs$sumAddSupps)),]
				}
				if(suppMethod=="minSum") {
					iqs <- iqs[which(iqs$sumAddSupps == min(iqs$sumAddSupps)),]
					iqs <- iqs[which(iqs$anzAddSupps == min(iqs$anzAddSupps)),]
				}
				if(suppMethod=="minSumLogs") {
					iqs <- iqs[which(iqs$sumAddSupps == min(log(1+iqs$sumAddSupps))),]
					iqs <- iqs[which(iqs$anzAddSupps == min(iqs$anzAddSupps)),]
				}
				
				# finally choose the suppression scheme
				if (randomResult == TRUE) {
					iqs <- iqs[sample(nrow(iqs), 1),]
				}
				else {
					iqs <- iqs[1,]
				}
				erg <- iqsInfo[[iqs$listInd]]
				erg$nrDiametral <- iqs$listInd
			}

			# we have a problem: there is no suppression scheme satisfying the nessessary interval protection
			else {
				cat("No suppression scheme satisfying the interval protection exist. -> The suppression scheme featuring the highest possible protection level is chosen!\n") 
				if(length(which(subtab$data$geh != "") == 0)) {
					#cat("Everything ok, the entire subtable is already suppressed!\n")
					erg <- NULL
				}
				else {
					iqs <- iqs[which(iqs$schutz == max(iqs$schutz)),]
					iqs <- iqs[which(iqs$anzAddSupps == min(iqs$anzAddSupps)),]
					iqs <- iqs[which(iqs$sumAddSupps == min(iqs$sumAddSupps)),]
					iqs <- iqs[1,]
					erg <- iqsInfo[[iqs$listInd]]
					erg$nrDiametral <- iqs$listInd
			   }
			}
			return(erg)
		}
		
		# suppress the values of the optimal suppression scheme
		suppressQuader <- function(subtab, optimalQuader) {
			t <- optimalQuader$quaderPosition[which(subtab$data$geh[optimalQuader$quaderPosition] == "")]
			# possibly a optimal scheme exists in which no additional values need to be suppressed, which, however, contains single cells
			if(length(t) > 0) {
				subtab$data[t,"geh"] <- "S"
			}
			return(subtab)
		}

		# find additional suppressions to protect single cells
		findAdditionalQuader <- function(info, optimalQuader) {
			iqs <- info$iqs
			indexList <- NULL

			quaderWithSingles <- do.call(rbind, lapply(iqs, '[', 'indikatorSingleItems'))
			quaderWithZeros <- do.call(rbind, lapply(iqs, '[', 'indikatorZeroItems'))

			zeros <- which(quaderWithZeros==FALSE)
			singles <- which(quaderWithSingles==FALSE)

			# normalized indices of single cells in the original optimal suppression scheme
			singleIndexOptQuader <- apply(optimalQuader$singleItems,1, function(x) { paste(x, collapse="") } )
			names(singleIndexOptQuader) <- NULL

			for (i in 1:length(iqs)) {
				# optimalQuader$nrDiametral is the number of the optimal index
				# this quader is not a possible additional suppression quader
				# if a scheme doesn not contain single cells, it is a potential additional suppression scheme
				if(i != optimalQuader$nrDiametral) {
					if(is.null(iqs[[i]]$singleItems)) {
						indexList <- cbind(indexList, i)
					}
					else {
						# the current suppression scheme contains single cells
						# check if these single values also exist in the oritinal suppression scheme
						# if not, it is a potential additional suppression scheme
						singleIndexAktQuader <- apply(iqs[[i]]$singleItems, 1, function(x) { paste(x, collapse="") } )
						names(singleIndexAktQuader) <- NULL

						if(length(which(singleIndexOptQuader %in% singleIndexAktQuader)) == 0) {
							indexList <- cbind(indexList, i)
						}
					}
				}
			}

			if(is.null(indexList)) {
				return(erg=NULL)
			}
			else {
				# check if there are potential suppression schemes that do not contain empty cells (they are preferred)
				erg <- list()
				ergZeros <- which(quaderWithZeros == FALSE)
				indexList <- as.vector(indexList)
				possibleAdditionalQuaderIndices <- which(ergZeros %in% indexList)
				# we create a new list with potential additional suppression schemes and chose an optimal quader afterwards
				if(length(possibleAdditionalQuaderIndices) > 0) {
					for(i in possibleAdditionalQuaderIndices) {
						erg[[i]] <- iqs[[i]]
					}
				}
				else {
					z <- 0
					for(i in indexList) {
						z <- z + 1
						erg[[z]] <- iqs[[i]]
					}
				}
				return(list(subtab=info$subtab, supp=info$supp, iqsInfo=erg))
			}
		}

		geh <- which(rownames(subtab$data) %in% subtab$supps2check)
		ind2 <- subtab$numberindexvars+2
		anzSupps <- anzSuppStart <- length(which(subtab$data$geh != ""))
		# we work out all values which needs to be protected sequentionally
		for (i in geh) {
			# g: index of the values which needs to be protected
			g <- as.numeric(as.character(subtab$data[i, (ind2+1):ncol(subtab$data)]))
			# supp: the potential diametral indices and information about the suppressed value 
			supp <- diametralIndex(subtab, g) 
			info <- calculateInformationForSuppValg(subtab, supp, allowZeros, protectionLevel)
			# only if the value is not yet protected
			if(!is.null(info)) {   
				optQuader <- findOptimalQuader(info, allowZeros, randomResult, suppMethod, protectionLevel)
				subtab <- suppressQuader(subtab, optQuader)
				# check single cells
				if(!is.null(optQuader)) {
					if(optQuader$indikatorSingleItems == TRUE) {
						cat("We need to find an additional suppression scheme.\n")
						addQ <- findAdditionalQuader(info, optQuader)
						if(!is.null(addQ)) {
							optQuaderAdditional <- findOptimalQuader(addQ, allowZeros=TRUE, randomResult=FALSE, suppMethod="minSupps", protectionLevel=0)
							if(!is.null(optQuaderAdditional)) {
								subtab <- suppressQuader(subtab, optQuaderAdditional)
							}
						}
						# no suppression scheme exists that does not contain single values -> we suppress the entire subtable
						else {
							subtab$data$geh[subtab$data$geh==""] <- "S"
						}
					}
				}
			}
		}
		subtab
	}

	# add necessary columns needed for GHMITER-approach
	fullData$data <- recodeIndexVars(fullData$data, fullData$indexvars)
	
	# if we check if secondary suppressions are backed up, we impose a very low protection level
	if (fullData$counter > 1) {
		protectionLevel <- 10
	}

	start <- as.numeric(Sys.time())
	SuppsStart <-  which(fullData$data$geh != "")
	fullData$data$geh[SuppsStart] <- "S"
	fullData$data$geh[fullData$supps2check] <- "P"	
	
	anzSuppStart <- length(SuppsStart)
	spl <- splitPrimarySupps(fullData)
	nc <- ncol(spl[[1]])
	# We run throuh all groups
	cat("The algorithm runs over all subgroups with primary/secondary suppressed values!\n")
	lSpl <- length(spl)
	for (i in 1:lSpl) {
		txtProgressBar(max = lSpl, initial = i, char = "=",style = 3)
		subtab <- calcSubset(fullData, spl[[i]])
		subtab <- algorithmGHMITER(subtab, allowZeros, randomResult, suppMethod, protectionLevel)
		ind <-  as.numeric(rownames(subset(subtab$data, subtab$data$geh=="S")))
		if(length(ind) > 0) {
			fullData$data[ind, "geh"] <- "S"
		}
	}
	cat("\n")
	
	SuppsEnd <- which(fullData$data$geh!= "")
	anzSuppEnd <- length(SuppsEnd)

	suppsNew <- SuppsEnd[-which(SuppsEnd %in% SuppsStart)]

	anzSuppSec <- anzSuppEnd - anzSuppStart
	end <- as.numeric(Sys.time())

	fullData$supps2check <- suppsNew

	return(list(fullData=fullData, anzSecSupp=anzSuppSec, time=end-start))
}
