\name{cellInfo}
\alias{cellInfo}
\title{cellInfo}
\description{cellInfo() calculates the anonymity state of a given cell.}
\usage{cellInfo(outObj, characteristics, varNames, return=FALSE) }
\arguments{
\item{outObj}{a data-object derived from primarySuppression (of class 'outObj')\cr
or protectTable (of class 'safeTable')}
\item{characteristics}{for each dimensional variable a given characteristic.\cr
It is important that the characteristics are in the same order as specified\cr
in input object 'varNames'}
\item{varNames}{vector of variable names of the dimensional variables which needs\cr
to be of exact same length as the number of dimensional variables and \cr
the length of input object 'characteristics'.}
\item{return}{if set to TRUE, an object containing the information calculated \cr
is returned, else the information is only printed}
}
\details{
if return is set to TRUE, a list-object is returned which containes interesting anonymity-information about the selected table cell.\cr
if return is set to FALSE, this information is only printed to the screen.
}
\value{manipulated data.}

\author{Bernhard Meindl}
\examples{
protectedData <- get(load(paste(searchpaths()[grep("sdcTable", searchpaths())], "/data/protectedData.RData", sep="")))

# check state of some cells...
a <- cellInfo(protectedData, c("011","B","0112","w"), c("V1", "V3", "V4", "V2"), return=TRUE) 
b <- cellInfo(protectedData, c("011","w","B","0112"), c("V1", "V2", "V3", "V4"), return=TRUE) 

print(a$problematicPrimaryCells)
print(b$problematicPrimaryCells)

# publishable cell
cellInfo(protectedData, c("A","013","0013","w"), c("V3", "V1", "V4", "V2"), return=FALSE) 
# a primary suppressed cell
cellInfo(protectedData, c("0013","m","A","021"), c("V4", "V2", "V3", "V1"), return=FALSE) 

# a secondary suppressed cell
cellInfo(protectedData, c("0011","m","A","011"), c("V4", "V2", "V3", "V1"), return=FALSE) 
}
\keyword{methods}