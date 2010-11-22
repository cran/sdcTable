protectTable <- function(fullData, method, ...){
	suppMethod = list(...)$suppMethod
	protectionLevel = list(...)$protectionLevel
	allowZeros = list(...)$allowZeros
	randomResult = list(...)$randomResult
	factorUp = list(...)$factorUp
	factorDown = list(...)$factorDown
	ub = list(...)$ub
	lb = list(...)$lb
	UPLPerc = list(...)$UPLPerc
	LPLPerc = list(...)$LPLPerc
	weight = list(...)$weight
	
	if(class(fullData) != "fullData")
		stop("Input object needs to be of class \"fullData\"\n")
	
	if(!method %in% c("HYPERCUBE", "HITAS")) 
		stop("Please choose a valid method! Choices are \"HYPERCUBE\" and \"HITAS\"!\n")

	# check if all necessary parameters are set
	if(method %in% c("HYPERCUBE")) {
		if(is.null(suppMethod))
			cat("the default value \"minSupps\" for parameter \"suppMethod\" will be used!\n")
		if(is.null(protectionLevel))
			cat("the default value of 80 for parameter \"protectionLevel\" will be used!\n")
		if(is.null(allowZeros))	
			cat("the default value \"TRUE\" for parameter \"allowZeros\" will be used!\n")
		if(is.null(randomResult))	
			cat("the default value \"FALSE\" for parameter \"randomResult\" will be used!\n")
	}

	if(method == "HITAS") {
		if(is.null(UPLPerc))
			cat("the default value of 15 for parameter \"UPLPerc\" will be used!\n")
		if(is.null(LPLPerc))
			cat("the default value of 15 for parameter \"LPLPerc\" will be used!\n")
		if(is.null(weight))
			cat("the default choice \"values\" for parameter \"weight\" will be used!\n")
	}

	indexvars <- fullData$indexvars
	
	erg <- list()
	erg$fullData <- fullData	
	
	supps <- NULL
	time <- NULL
	
	counter <- 1
	
	# original primary suppressed cells
	origPrimarySuppressions <- fullData$supps2check

	if(method=="HITAS") {
		time <- Sys.time()
		erg <- processTableHITAS(erg$fullData, ...) 
		time <- Sys.time() - time
	
		erg$time <- time
		erg$counter <- 1
		erg$method <- method
		erg$totSupps <- length(which(erg$fullData$data$geh=="S"))
	}
	else {  
		ind <- FALSE	
	    while(ind==FALSE) {					
	        cat("Cycle",counter,"to protect values is now started ...\n")
			if(method=="HYPERCUBE"){
				erg$fullData$counter <- counter
				erg <- processTableHYPERCUBE(erg$fullData, ...)
			}
			
			time <- cbind(time, as.numeric(erg$time))
			supps <- cbind(supps, erg$anzSecSupp)
	        if(erg$anzSecSupp == 0)
	            ind <- TRUE
			counter <- counter + 1
	    }	    
	    
		# which cells are primary/secondary suppressed?
		suppressions <- which(erg$fullData$data$geh != "")
		secondarySuppressions <- suppressions[which(! suppressions %in% origPrimarySuppressions)]
		
		erg$fullData$data$geh[origPrimarySuppressions] <- "P"
		erg$fullData$data$geh[secondarySuppressions] <- "S"	
		
		## NEW: additionally suppress all dimensions with only one sub-level
		for(i in 1:length(indexvars)) {
			dims <- erg$fullData$dimensions[[i]]
			for(j in 1:length(dims)) {
				if(length(dims[[j]])==2) {					
					indTest <- which(erg$fullData$data[,indexvars[i]]%in% dims[[j]])
					ss <- erg$fullData$data[indTest,]
					if(length(indexvars) > 1) {
						spl <- split(ss, ss[,indexvars[-i]])		
						for(z in 1:length(spl))
							if(length(unique(spl[[z]][,"geh"]))==2)
								erg$fullData$data[rownames(spl[[z]]),"geh"] <- sort(spl[[z]][,"geh"])[2]
					}
					else {
						if(length(unique(ss$geh))==2)
							erg$fullData$data[rownames(ss),"geh"] <- sort(ss[,"geh"])[2]
					}
				}					
			}
		}
		
		
		erg$counter <- counter
	    erg$time <- sum(time)
		erg$method <- method
		erg$totSupps <- sum(supps)	
	}	
	class(erg) <- "safeTable"
    return(erg)
}