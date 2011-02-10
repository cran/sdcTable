\name{primarySuppression}
\alias{primarySuppression}
\title{primarySuppression}
\description{primarySuppression() allows to modify an object created by createFullTable() to specfiy the cells that need to be protected by one of three popular primary suppression rules.}
\usage{primarySuppression(outObj, suppRule_Freq=NULL, suppRule_P=NULL, suppRule_NK=NULL)}
\arguments{
\item{outObj}{an object created with createFullTable().}
\item{suppRule_Freq}{a vector of length 2. the first element specifies the parameter 'n' of the n-rule meaning that all cells having less than n elements. if the second element of suppRule_Freq equals 0, then empty cells will not be considered primary unsafe. }
\item{suppRule_P}{vector of length 1. the first element specifices parameter 'p' of the (p)-percent suppression rule.}
\item{suppRule_NK}{vector of length 2. the first element specifices parameter 'n', the second element specifies parameter 'k' of the (n,k)-suppression rule.}
}
\value{manipulated data.}

\author{Bernhard Meindl}
\examples{
datObj <- get(load(paste(searchpaths()[grep("sdcTable", searchpaths())], "/data/testObj.RData", sep="")))

# all cells with less than 3 units should be primary unsafe
suppRule_Freq <- c(3, 0)
datObj1 <- primarySuppression(datObj, suppRule_Freq=suppRule_Freq)

# p-percent rule with p=80
suppRule_P 	<- 80
datObj1 <- primarySuppression(datObj, suppRule_P=suppRule_P)
}
\keyword{methods}