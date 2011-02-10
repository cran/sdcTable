prepareInput <- function(dat, filenames=NULL, hierFrames=NULL, freqVar=NULL, numVar=NULL, weightVar=NULL, sampWeightVar=NULL, suppRule_Freq=NULL, suppRule_P=NULL, suppRule_NK=NULL) {
	# use filename without extension as variable-name. needs to match
	# TODO: generalize
	if ( is.null(filenames) & is.null(hierFrames) )
		stop("please specify either filenames or hierFrames!\n")
	if ( !is.null(filenames) & !is.null(hierFrames) )
		stop("please specify filenames or hierFrames!\n")
	if ( !is.null(filenames) ) 
		varnames <- as.vector(sapply(filenames, function(x) { strsplit(tail(unlist(strsplit(x, "/")),1),"\\.")[[1]] } )[1,])
	if ( !is.null(hierFrames) ) {
		if ( !is.list(hierFrames) )
			stop("hierFrames needs to be a list of data.frames with list-names being variable-names existing in the input data-set!\n")
		else
			varnames <- names(hierFrames)
	}
	
	if ( is.data.frame(dat) ) 
		vnames <- colnames(dat)
	else if ( is.list(dat) ) {
		vnames <- names(dat)
		dat <- as.data.frame(dat)
	}
	else {
		stop("input object data needs to be a data.frame or list!\n")
	}
	
	# 0: calc position
	varPos <- match(varnames, vnames)
		
	# 1: calc levelObj
	levelObj <- list()
	for ( i in 1:length(varPos) ) {
		if ( !is.null(filenames) )
			levelObj[[i]] <- calcDimInfos(dat, file=filenames[i], dataframe=NULL, vName=varnames[i]) 
		else
			levelObj[[i]] <- calcDimInfos(dat, file=NULL, dataframe=hierFrames[[i]], vName=varnames[i]) 
	}	
	
	# 2: calcFullTable with respect to the case of microdata or aggregated data
	outObj <- calcFullTable(dat, levelObj, freqVar, numVar, weightVar, sampWeightVar)	
		
	# 3: are all levels (minimal=TRUE) available? 
	# sollte in calcFullTable schon erledigt sein!
		
	# 4: optionally: do primary suppression	
	if ( any(!is.null(c(suppRule_Freq, suppRule_P, suppRule_NK))) )
		outObj <- primarySuppression(outObj, suppRule_Freq, suppRule_P, suppRule_NK)
	
	# 5) TODO: optionally possibility to change bounds using setBounds()	
	outObj	
}
