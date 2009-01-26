\name{summary.safeTable}
\alias{summary.safeTable}
\title{summary.safeTable}
\description{summary method for objects of class safeTable.}
\usage{\method{summary}{safeTable}(object, ...)}
\arguments{
	\item{object}{object from class safeTable}
	\item{...}{additional parameters. Not used yet.}
}
\details{object is an object of class safeTable. The summary functions returns several statistics from the anonymisation process.}
\value{Manipulated data.}
\author{Bernhard Meindl}
\examples{
	\dontrun{
		data(exampleFullData)
		erg <- protectTable(exampleFullData, method="OPT")
		class(erg)
		summary(erg)
	}
}
\keyword{print}