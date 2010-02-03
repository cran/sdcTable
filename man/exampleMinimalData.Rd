\name{exampleMinimalData}
\alias{exampleMinimalData}
\docType{data}
\title{minimal data set used as input to generate complete data set.}
\description{a data set containing a total of 3 dimensional variables 
and the corresponding cell values. Only those characteristics of the 
dimensional variables are available in exampleMinimalData that cannot be 
calculated as linear combinations from other characteristics (sums).}
\usage{data(exampleMinimalData)}
\format{
A data.frame containing a column for each dimensional variable and a column for the corresponding cell value.
\describe{
\item{\code{Var1}}{characteristic of dimensional variable 1.}
\item{\code{Var2}}{characteristic of dimensional variable 2.}
\item{\code{Var3}}{characteristic of dimensional variable 3.}
\item{\code{val}}{cell value of cell given by characteristics of Var1, Var2 and Var3.}
}
}
\examples{
\dontrun{
data(exampleMinimalData)
head(exampleMinimalData)
dim(exampleMinimalData)
}
}
\keyword{datasets}