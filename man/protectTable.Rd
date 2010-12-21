\name{protectTable}
\alias{protectTable}
\title{protectTable}
\description{Protect hierarchical tabular data using secondary suppression.}
\usage{protectTable (fullData, method, ...)}
\arguments{
  \item{fullData}{object from class fullData.}
  \item{method}{either HYPERCUBE, HITAS}
  \item{...}{additional input parameters.}
}
\details{
The additional input parameters depend on secondary suppression algorithm
(HITAS or HYPERCUBE are currently supported). \cr
If one chooses "HYPERCUBE", it is possible to specify the following parameters:\cr
a) allowZeros: TRUE if cells with frequency of zero should be allowed to be part of the suppression scheme, FALSE otherwise \cr
b) randomResult: if TRUE the procedure will return different protection schemes in
case of several possibilities (which are identical with respect to the loss-function)\cr
c) suppMethod: the algorithm minimizes the number of additional suppressions ("minSupps"),
the sum of additional suppressed values ("minSum") or the logarithmized sum of 
additional suppressed cells ("minSumLogs")\cr
d) protectionLevel: specifies the desired protection level in percent\cr\cr
If "HITAS" is chosen, the following parameters can be set:\cr
a) ub: upper bounds for cell values assumed to be known by attackers. If not specified, cell value * 100 is assumed.\cr
b) lb: lower bounds for cell values assumed to be known by attackers. If not specified, 0 is assumed for each cell (non-negative table).\cr
c) UPLPerc: percentage of (upper) protection level for each primary sensitive cell.\cr
d) LPLPerc: percentage of (lower) protection level for each primary sensitive cell.\cr
e) weight: parameter to use for objective function. Choices are "values" and "logs".
}
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
minDat <- expand.grid(V1=V1, V2=V2)
minDat$value <- rpois(nrow(minDat), 7)
minimalDat <- minDat
l <- list()
l[[1]] <- l[[2]] <- c(1,1)
tableVars=c("V1","V2")
freqRule <- c(3, 0) # primary suppress values >0 & <= 3
numVar <- NULL
freqVar="value"
# creating the full data set needed for protectTable()
fullDat <- createFullData(
	minimalDat=minDat, 
	tableVars=tableVars, 
	numVar=numVar, 
	freqVar=freqVar,
	suppRule_Freq=freqRule,
	l=l) 
	
# protecting the data using HYPERCUBE and HITAS approach
protData1 <- protectTable(fullDat, method="HYPERCUBE")
protData2 <- protectTable(fullDat, method="HITAS")

summary(protData1); summary(protData2)

### example with existing data-set
data(exampleFullData)
result1 <- protectTable(
			exampleFullData, 
			method="HYPERCUBE", 
			allowZeros=FALSE, 
			suppMethod="minSum")
result2 <- protectTable(
			exampleFullData, 
			method="HITAS", 
			LPLPerc=15, 
			UPLPerc=20)
summary(result1)
summary(result2)
}	
}
\keyword{methods}