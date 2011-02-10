\name{checkSuppressionPattern}
\alias{checkSuppressionPattern}
\title{checkSuppressionPattern}
\description{checkSuppressionPattern() calculates the anonymity state of a given cell.}
\usage{checkSuppressionPattern(outObj, pattern, debug=TRUE, all=FALSE, stop=TRUE)}
\arguments{
\item{outObj}{a data-object of class 'safeTable' derived from protectTable()}
\item{pattern}{a vector consisting of TRUE|FALSE only with the same length\cr
as the total number of possible table cells in the same order as in outObj\cr
specifying for each table cell if the cell is part of the suppression scheme (TRUE)
or not (FALSE).}
\item{debug}{if debug)==TRUE, for each checked cell the calculated interval is printed.\cr
If debug==FALSE, nothing is printed.}
\item{all}{if FALSE, only primary suppressed cells are checked. If all==TRUE, then \cr
also secondary suppressed cells are checked if they cannot be recalculated.}
\item{stop}{if TRUE, the checking procedure stops as soon as at least one cell is not safe.}
}
\value{returns a list with three objects. 'validPattern' is 'TRUE' if 'pattern' is\cr
a valid suppression pattern for the given input data and 'FALSE' otherwise. 'limits' is a \cr
list containing to calculated lower and upper bounds. The third object 'indices' returns the\cr
indices of the cells that have been checked.}

\author{Bernhard Meindl}
\examples{
\dontrun{ 
protectedData <- get(load(paste(searchpaths()[grep("sdcTable", searchpaths())], "/data/protectedData.RData", sep="")))
suppVec <- protectedData$outObj$status

pattern <- rep(FALSE, length(suppVec))
pattern[!is.na(match(suppVec, c("u","x")))] <- TRUE

c1 <- checkSuppressionPattern(protectedData, pattern, debug=FALSE, all=FALSE)

# remove a cell from the suppression pattern
pattern[which(suppVec=="x")[1]] <- FALSE
c2 <- checkSuppressionPattern(protectedData, pattern, debug=FALSE, all=FALSE, stop=FALSE)

print(c1$validPattern)
print(c2$validPattern)
}
}
\keyword{methods}