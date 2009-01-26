\name{exampleFullData}
\alias{exampleFullData}
\docType{data}
\title{Example complete data set used as input for various methods}
\description{exampleFullData is a data set of class 'fullData'}
\usage{data(exampleFullData)}
\format{
	A list containing several necessary information.
	\describe{
    	\item{\code{data}}{data set containing a column for each dimension, a column representing the numerical value of each cell and a column which indicates if a cell is marked as sensitive (P) or not.}
    	\item{\code{dims}}{a list representing the structure of each dimensional variable.}
    	\item{\code{indexvars}}{a vector representing the columns of the dimensional variables in data set.}
    	\item{\code{supps2check}}{a numeric vector of indices representing cells that need to be protected.}
		\item{\code{numberindexvars}}{the total number of dimensional variables.}
	}
}
\examples{
	\dontrun{
		data(exampleFullData)
		head(exampleFullData$data)
		dim(exampleFullData$data)
	}
}
\keyword{datasets}