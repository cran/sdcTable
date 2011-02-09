\name{protectLinkedTables}
\alias{protectLinkedTables}
\title{protectLinkedTables}
\description{protectLinkedTables() allows to protect linked data-object. 'Linked' means that e.g tables feature at least one common cell. Therefore, if the aim is to protect the data under consideration, it is neccessary to take special care of the common cells since these cells need to have the same status (suppressed or not suppressed) after the protection procedure.\cr
The common cells are specified using the input object 'commonCells' which needs to be a list in a specific format.\cr
The algorithm iteratively protects the data-sets and checks if the stop criterion (all common cells have the same suppression status) is fulfilled. If so, the procedure stops. If at least one common cell has different status, the this cell is set to primary suppressed in the other dataset and the protection procedure starts again. Please note, that this iterative algorithm may lead to significant over-suppression. }
\usage{protectLinkedTables(inputObj1, inputObj2, commonCells, method="HITAS", weight=NULL)}
\arguments{
\item{ inputObj1 }{ a data-object created by calcFullTable() and primarySuppression() }
\item{ inputObj2 }{ a data-object created by calcFullTable() and primarySuppression() }
\item{ commonCells }{ a list object specifying 'common cells' between inputObj1 and inputObj2.\cr
Each list element of 'commonCells' has to be a list element too. For each of these list-elements there are two possible choices.\cr
The first choice has to be used if a dimension exists in both input-data objects. In this case a list-object of length 3 has to be specified. The first element specifies the position of the variable under consideration in the first dataset, the second list element its position in the second (micro)dataset (see example below). The third element consists of the keyword 'ALL' which tells protectLinkedTables that the variable under consideration is equal in both input datasets.\cr
The second possibility has to be used if different variables have common (sub)totals or cells. In this case a list-object of length 4 has to be specified. Element 1 and 2 specify the position of the variables in data-set 1 and two (just as described before). The third list-element is a vector of the characteristics of the variable in dataset 1 for which identical entries exist in dataset 2. The fourth list-element is a vector too, specifying the corresponding characteristics for the dimensional variable from dataset 2 (please have a look at the example provided below!). }
\item{ method }{ choice of suppression algorithm. Currently 'HITAS' and 'HYPERCUBE' are valid choices. }
\item{ weight }{ currently not used. }
}
\value{ manipulated data. }
\author{Bernhard Meindl}
\examples{
# generate some micro-data
# NOTE: we do this in a way that EcoOld and EcoNew have common cells when
# aggregating over the other dimensions.
N <- 200
Region <- sample(c("01","02"), N, replace=TRUE)
Sex <- sample(c("01","02"), N, replace=TRUE)
EcoOld <- sample(c("011","012","021","022"), N, replace=TRUE)
microDat <- data.frame(Region,Sex,EcoOld, EcoNew=NA)

spl <- split(microDat, apply(microDat[,1:2], 1, paste, collapse=""))
for ( i in 1:length(spl) ) {
	ind1 <- which(substr(spl[[i]]$NaceOld,1,2)=="01")		
	ind2 <- setdiff(1:nrow(spl[[i]]), ind1)
	if ( length(ind1) > 0 )
		spl[[i]]$EcoNew[ind1] <- sample(c("011", "012","013"), length(ind1), replace=TRUE)
	if ( length(ind2) > 0 )
		spl[[i]]$EcoNew[ind2] <- sample(c("021","022","023"), length(ind2), replace=TRUE)
}
microDat <- do.call("rbind", spl)
rownames(microDat) <- 1:N
microDat$numVal <- abs(round(rnorm(N, 500, 200),2))
	
microDat1 <- microDat[,c(2,3,5)] 	# Sex, EcoOld and numVal
microDat2 <- microDat[,c(1,2,4,5)]	# Region, Sex, EcoNew and newVal
		
# Region: exists only in microDat2
df1 <- data.frame(h=c("@@","@@"), l=c("R1","R2"))
dim1b <- calcDimInfos(microDat2, file=NULL, dataframe=df1, vName="Region") 
	
# Sex: exists in microDat1 and microDat2
df2 <- data.frame(h=c("@@","@@"), l=c("m","f"))
dim2a <- calcDimInfos(microDat1, file=NULL, dataframe=df2, vName="Sex") 
dim2b <- calcDimInfos(microDat2, file=NULL, dataframe=df2, vName="Sex") 
	
# Economic classification: (old version, exists only in microDat1)
df31 <- data.frame(
		h=c("@@","@@@","@@@","@@","@@@","@@@"), 
		l=c("A","Aa","Ab","B","Ba","Bb"))
dim31a <- calcDimInfos(microDat1, file=NULL, dataframe=df31, vName="EcoOld") 

#Economic classification: (new version, exists only in microDat2)
df32 <- data.frame(
		h=c("@@","@@@","@@@","@@@","@@","@@@","@@@","@@@"), 
		l=c("C","Ca","Cb","Cc","D","Da","Db","Dc"))
dim32b <- calcDimInfos(microDat2, file=NULL, dataframe=df32, vName="EcoNew") 
		
# the complete levelObjects
levelObj1 <- list(dim2a, dim31a)			# Sex, EcoOld	
levelObj2 <- list(dim1b, dim2b, dim32b)		# Region, Sex, EcoNew

numVar <- "numVal"	# the variable name of the numeric variable
suppRule_Freq <- c(15, 0)	# a simple rule for primary suppression

inputObj1 <- calcFullTable(microDat1, levelObj1, numVar)
inputObj1 <- primarySuppression(inputObj1, suppRule_Freq=suppRule_Freq)

inputObj2 <- calcFullTable(microDat2, levelObj2, numVar)
inputObj2 <- primarySuppression(inputObj2, suppRule_Freq=suppRule_Freq)

# specifiying common cells
commonCells <- list()
# variable "Sex"
commonCells[[1]] <- list()
commonCells[[1]][[1]] <- 1 # first column in microDat1
commonCells[[1]][[2]] <- 2 # second column in microDat2
commonCells[[1]][[3]] <- "ALL" # Sex has equal characteristics on both datasets
# Economic classification
commonCells[[2]] <- list()
commonCells[[2]][[1]] <- 2	# economic classification (old version) is second column in microDat1
commonCells[[2]][[2]] <- 3	# economic classification (new version) is third column in microDat2
commonCells[[2]][[3]] <- c("A","B")	# vector of common characteristics: A and B in ecoOld
commonCells[[2]][[4]] <- c("C","D") # correspond to C and D in ecoNew!
out <- protectLinkedTables(inputObj1, inputObj2, commonCells, method="HITAS")
print(str(out))
}
\keyword{methods}