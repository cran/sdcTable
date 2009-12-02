\name{protectTable}
\alias{protectTable}
\title{protectTable}
\description{Wrapper-function to call several algorithms for secondary cell suppression of tabular data.}
\usage{protectTable (fullData, method, ...)}
\arguments{
  \item{fullData}{object from class fullData.}
  \item{method}{either HYPERCUBE, HITAS}
  \item{...}{additional input parameters.}
}
\details{For correct input parameters for function protectTable() please have a look at Please have a look at the corresponding functions processTableHYPERCUBE() or rather processTableHITAS().}
\value{Manipulated data.}
\references{
	Repsilber, D. (1999). Das Quaderverfahren. In: Forum der Bundesstatistik, Band 31/1999. \\
	de Wolf, P.P (2002). HiTaS: A Heuristic Approach to Cell Suppression in Hierarchical Tables. In: Domingo-Ferrer, J. (Hrsg.): Inference Control in Statistical Databases. Vol. 2316. \\
	Fischetti, M., Salazar, J.J. (2000). Models and Algorithms for Optimizing Cell Suppression in Tabular Data with Linear Constraints. In: Journal of the American Statistical Association 95, 916-928.
}
\author{Bernhard Meindl}
\examples{
	\dontrun{
		### very simple example: protect a 2-dimensional table with only one marginal sum for each dimensional variable
		V1 <- c("01", "02", "03")
		V2 <- c("01", "02", "03","04", "05")
		minDat <- expand.grid(V1, V2)
		minDat$value <- rpois(nrow(minDat), 7)
		
		# creating the full data set needed for protectTable()
		fullDat <- createFullData (minDat, indexvars=1:2, l=list(c(1,1), c(1,1)), suppVals=TRUE, suppLimit=3, suppZeros=FALSE)
		
		# protecting the data using HYPERCUBE and HITAS approach
		protData1 <- protectTable(fullDat, method="HYPERCUBE")
		protData2 <- protectTable(fullDat, method="HITAS")
		
		summary(protData1); summary(protData2)
		
		### example with existing data-set
		data(exampleFullData)
		result1 <- protectTable(exampleFullData, method="HYPERCUBE", allowZeros=FALSE, suppMethod="minSum")
		result2 <- protectTable(exampleFullData, method="HITAS", LPLPerc=15, UPLPerc=20)
		summary(result1)
		summary(result2)
	}	
}
\keyword{methods}