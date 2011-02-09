\name{calcFullTable}
\alias{calcFullTable}
\title{calcFullTable}
\description{calcFullTable() takes input data (both microdata and already aggregated data are possible together with information about all dimensional variables (defined by 'levelObj') and returns a list-object that contains all the information that is neccessary to protect the input data. calcFullTable() calculates the complete hierarchical structure of the input-data, sets default bounds for the required cell-protection levels (which might be changed with setBounds(). If the parameter 'numVar' is specified, the corresponding column is treated as numerical variable (and not as the number of units contributing to a cell). }
\usage{calcFullTable(dataset, levelObj, freqVar=NULL, numVar=NULL, weightVar=NULL, sampWeight=NULL)}
\arguments{
\item{dataset}{ a dataframe or matrix. The input data consists of columns defining the dimensions and (optionally) a column with a numeric variable (defined by 'numVar'). Furthermore (if the input data are already aggregated data) the user may specify a column showing the number of units contributing to a given cell. }
\item{levelObj}{ a list of objects (for each dimensional-variable). Each list-element has to be created with function calcDimInfos(). }
\item{freqVar}{ if not NULL, 'freqVar' specifies the column name of the count variable in 'dataset'. If 'NULL' it is assumed that we are dealing with microdata and therefore each row in dataset is assigned frequency 1.}
\item{numVar}{ if not NULL, 'numVar' specifies the column name of the numerical variable in 'dataset'. }
\item{weightVar}{ if given, the variable name of a variable containing weights used in the CSP master-problem.}
\item{sampWeight}{ if given, the variable name of a variable used for weighting purposes.}
}
\value{ manipulated data. }
\author{Bernhard Meindl}
\examples{
	# the level object (information about dimensions)
	levelObj <- get(load(paste(searchpaths()[grep("sdcTable", searchpaths())], "/data/levelObj.RData", sep="")))
	
	# micro data
	microDat <- get(load(paste(searchpaths()[grep("sdcTable", searchpaths())], "/data/microDat.RData", sep="")))
	print(head(microDat))
	
	# aggregated data
	aggDat <- get(load(paste(searchpaths()[grep("sdcTable", searchpaths())], "/data/aggregatedDat.RData", sep="")))
	print(head(aggDat))	
	
	out1 <- calcFullTable(microDat, levelObj, freqVar=NULL, numVar="numVal", weightVar=NULL)
	out2 <- calcFullTable(aggDat, levelObj, freqVar="Freq", numVar="numVal", weightVar=NULL)
	
	# compare
	print(identical(out1$fullTabObj, out2$fullTabObj))
	print(str(out1))
	print(str(out2))
	
	aggDat$sampWeight <- sample(rpois(nrow(aggDat), 10), replace=TRUE)
	out3 <- calcFullTable(aggDat, levelObj, freqVar="Freq", numVar="numVal", sampWeight="sampWeight")
}
\keyword{methods}