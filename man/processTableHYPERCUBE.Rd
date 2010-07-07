\name{processTableHYPERCUBE}
\alias{processTableHYPERCUBE}
\title{processTableHYPERCUBE}
\description{HYPERCUBE - algorithm for secondary cell suppression.}
\usage{processTableHYPERCUBE (fullData, allowZeros=TRUE, randomResult=FALSE, suppMethod="minSupps", protectionLevel=80)}
\arguments{
	\item{fullData}{object from class fullData.}
	\item{allowZeros}{should empty cells be allowed to be part of the suppression scheme?}
	\item{randomResult}{if several possible suppression schemes exist, always use the same or pick a random one?}
	\item{suppMethod}{should the algorithm minimize the number of suppressions (minSupps), minimize the suppressed sum of values (minSum) or the sum of logarithmized values (minSumLogs)?}
	\item{protectionLevel}{protection level in percent.}
}
\details{Have a look at the link given below.}
\value{Manipulated data.}
\references{Repsilber, D. (1999). Das Quaderverfahren. In: Forum der Bundesstatistik, Band 31/1999}
\author{Bernhard Meindl}
\note{The algorithm runs only once through the entire data set. To make completely sure that also secondary suppressions are enough protected one should use function protectTable(x, method="GHMITER") which loops through the data set to ensure protection.}
\examples{
\dontrun{
data(exampleFullData)
erg <- processTableHYPERCUBE(exampleFullData, suppMethod="minSUM")
summary(erg)
}
}
\keyword{methods}

