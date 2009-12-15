\name{createFullData}
\alias{createFullData}
\title{createFullData}
\description{Function to generate an object of class fullData \cr 
			 which can then be used as input object for several \cr
			 other functions in sdcTable such as protectTable.}
\usage{createFullData(minimalDat=NULL, microDat=NULL, tableVars=NULL, numVar=NULL, suppRule_Freq=NULL, suppRule_P=NULL, suppRule_NK=NULL, l=l)}
\arguments{
  \item{ minimalDat }{a data.frame containing a column for each hierarchical \cr variable in a specific coding and a column specifying either the frequency of the cell or some aggregated numerical value.}
  \item{ microDat }{a data.frame containing a column for each hierarchical variable in a specific coding and (optionally) a column for the corresponding numerical values.}
  \item{ tableVars }{a vector containing the variable names of the dimensional variables in either microDat or minimalDat.}
  \item{ numVar }{either NULL or a variable name of the numeric variable of the table.}
  \item{ suppRule_Freq }{NULL or a vector of length 2. If specified, the first vector element gives the threshold at which cells should be suppressed while the second (TRUE|FALSE) tells, if zero-cells should be suppressed or not.}
  \item{ suppRule_P }{either NULL or - if specified - a vector of length 1 specifying the p-parameter used in the p\%-rule. This parameter is only evaluated, if microdata have been provided as input.}
  \item{ suppRule_NK }{either NULL or - if specified - a vector of length 2. The first element	specifying the n-parameter, the second element specifiying the k-parameter used in the nk-dominance-rule. This parameter is only evaluated, if microdata have been provided as input.}  
  \item{ l }{a list which element i defining the hierarchical structure of the i-th hierarchical variable.}
}
\details{Have a look at the link given below.}
\value{object of class fullData}
\author{Bernhard Meindl}
\note{fix me: LOTS! eg.,simplify the generation of object of class fullData? }
\examples{
	\dontrun{
	# generate micro-data
	genMicroData <- function(N) {
		V1 <- sample(c("011", "012","013","021","022"), N, replace=T)
		V2 <- sample(c("01", "02"), N, replace=T)
		V3 <- sample(c("01", "02"), N, replace=T)
		V4 <- sample(c("0011", "0012", "0013","0021","0031","0041","0042", "0051","0061","0062","0071","0072","0081","0091","0101","0102","0111","0112"), N, replace=T)
		
		microDat <- data.frame(V1=V1,V2=V2,V3=V3,V4=V4)
		microDat$numVal <- abs(round(rnorm(N, 500, 200),2))
		sInd <- sample(floor(N/20))
		microDat$numVal[sInd] <- abs(round(rnorm(sInd, 100000, 200),2))
		microDat
	}	
	# generate minimal data set (a table without subtotals)
	genMinimalDat <- function(micro, tableVars, numVar) {
		minimalDat <- as.data.frame(table(microDat[,which(colnames(microDat) \%in\% tableVars)]))
		indexvars <- which(colnames(micro) \%in\% tableVars)
		if(!is.null(numVar)) {
			if(length(tableVars)>1) {
				minimalDat$fac <- apply(minimalDat[,indexvars], 1, function(x) { paste(x, collapse="") } )
				fac <- apply(microDat[,indexvars], 1, function(x) { paste(x, collapse="") } )
			}
			else {
				minimalDat$fac <- minimalDat[,indexvars]
				fac <- microDat[,indexvars]		
			}
		}
		agg <- aggregate(microDat[,which(colnames(microDat)==numVar)], by=list(fac), sum)
		colnames(agg) <- c("fac", numVar)
		minimalDat <- merge(minimalDat, agg, all.x=T)
		minimalDat <- minimalDat[,-which(colnames(minimalDat)=="fac")]
		colnames(minimalDat)[1:length(tableVars)] <- tableVars
		minimalDat
	}		
		
	microDat <- genMicroData(2000)
	tableVars=c("V1","V2","V3","V4")
	numVar <- "numVal"
	minimalDat <- genMinimalDat(microDat, tableVars, numVar)
	
	suppRule_Freq <- c(3, FALSE)
	suppRule_P <- c(80) 
	#suppRule_NK <- c(2,75)
	
	l <- list()
	l[[1]] <- c(1,1,1)
	l[[2]] <- c(1,1)
	l[[3]] <- c(1,1)
	l[[4]] <- c(1,2,1)	
	
	# minimal Data
	# note: specifying suppRule_P is useless!
	result1 <- createFullData(minimalDat=minimalDat, tableVars=tableVars, numVar=numVar, suppRule_Freq=suppRule_Freq, suppRule_P=suppRule_P, l=l)
	print(class(result1))
	print(head(result1$data))
	
	# micro Data
	result2 <- createFullData(microDat=microDat, tableVars=tableVars, numVar=numVar, suppRule_Freq=suppRule_Freq, suppRule_P=suppRule_P, l=l)
	print(class(result2))
	print(head(result2$data))

	# result1 and result2 can now be used as input parameters for protectTable()
	}
}
\keyword{methods}