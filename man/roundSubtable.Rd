\name{roundSubtable}
\alias{roundSubtable}
\title{roundSubtable}
\description{Wrapper-function to call several algorithms for rounding algorithms.}
\usage{roundSubtable(subtab, indexvars=NULL, base=5, maxS=NULL, method=NULL)}
\arguments{
  \item{subtab}{subtable derived from object of class fullData.}
  \item{indexvars}{position of dimensional variables (needed only for controlled rounding procedure!).}
  \item{base}{the base-number to which cell values should be rounded.}
  \item{maxS}{maximum number of multiples of the base numbers at which a cell value is allowed to change (rounded up or down).}
  \item{method}{either \"simple\", \"random\" or \"controlled\" are valid choices.}
}
\details{Depending on the choice for parameter method, simple (or rather conventional), random or controlled rounding will be done.}
\value{Manipulated data.}
\references{
	Salazar, J.J., Bycroft C., Staggemeier A. (2005). Controlled Rounding Implementation. Supporting paper for Joint UNECE/EUROSTAT work session on statistical data confidentiality.
}
\author{Bernhard Meindl}
\examples{
\dontrun{
V1 <- c(rep("01", 4), rep("02", 4), rep("03", 4), rep("00", 4))
V2 <- rep(c("01", "02", "03", "00" ), 4)
vals <- c(20, 50, 10, 80, 8, 19, 22, 49, 17, 32, 12, 61, 45, 101, 44, 190)
subtab <- data.frame(V1, V2, vals)
indexvars <- 1:2
		
sr <- roundSubtable(subtab, base=3, method="simple")
rr <- roundSubtable(subtab, base=4, method="random")
cr <- roundSubtable(subtab, indexvars, base=5, maxS=3, method="controlled")
print(cbind(subtab, simpleRounding=sr$val, randomRounding=rr$val, controlledRounding=cr$val))
}	
}
\keyword{methods}