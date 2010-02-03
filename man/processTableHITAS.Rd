\name{processTableHITAS}
\alias{processTableHITAS}
\title{processTableHITAS}
\description{HITAS - algorithm for secondary cell suppression using optimal 
suppression algorithm.}
\usage{processTableHITAS (fullData, ub=NULL, lb=NULL, UPLPerc=35, LPLPerc=25, weight="values")}
\arguments{
	\item{fullData}{object from class fullData.}
  	\item{ub}{probably known upper bounds for cell values. If not specified, cellvalue * 100 is assumed. }
  	\item{lb}{probably known lower bounds for cell values. If not specified, 0 is assumed for each cell (non-negative table). }
  	\item{UPLPerc}{percentage of (upper) protection for each cell. Maximum possible value calculated for any sensitive cell must be greater or equal the cell value plus ULPPerc \%.}
  	\item{LPLPerc}{percentage of (lower) protection for each cell. Minimum possible value calculated for any sensitive cell must be less or equal the cell value minus LPLPerc \%.}
	\item{weight}{parameter to use for objective function. Possible choices are values and logs.}
}
\details{Have a look at the links given below.}
\value{Manipulated data.}
\references{	
	de Wolf, P.P (2002). HiTaS: A Heuristic Approach to Cell Suppression in Hierarchical Tables. In: Domingo-Ferrer, J. (Hrsg.): Inference Control in Statistical Databases. Vol. 2316. \\
	Fischetti, M., Salazar, J.J. (2000). Models and Algorithms for Optimizing Cell Suppression in Tabular Data with Linear Constraints. In: Journal of the American Statistical Association 95, 916-928.
}
\author{Bernhard Meindl}
\note{processTableHITAS() protects hierarchical tabular data using the HiTaS approach.}
\examples{
\dontrun{
data(exampleFullData)
erg <- processTableHITAS(exampleFullData, UPLPerc=15, LPLPerc=15)
summary(erg)
}
}
\keyword{methods}