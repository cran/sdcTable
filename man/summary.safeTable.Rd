\name{summary.safeTable}
\alias{summary.safeTable}
\title{summary.safeTable}
\description{summary method for objects of class safeTable.}
\usage{\method{summary}{safeTable}(object, ...)}
\arguments{
\item{object}{object from class safeTable}
\item{...}{additional parameters. Not used yet.}
}
\details{object is an object of class safeTable. The summary functions 
returns several statistics from the anonymisation process.}
\value{Manipulated data.}
\author{Bernhard Meindl}
\examples{
N <- 100
# generate micro-data
V1 <- sample(c("011", "012","013","021","022","023","024"), N, replace=TRUE)
V2 <- sample(c("01", "02"), N, replace=TRUE)
V3 <- sample(c("01", "02"), N, replace=TRUE)

microDat <- data.frame(V1=V1,V2=V2,V3=V3, numVal=abs(round(rnorm(N, 500, 200),2)))

# dimensional information (level1-level4
h1 <- c("@@", "@@@","@@@","@@@","@@", "@@@","@@@","@@@","@@@")
l1 <- c("010", "011", "012","013","020", "021","022","023","024")
df1 <- data.frame(h=h1, l=l1) #V1

h2 <- c("@@", "@@")
l2 <- c("m", "w")
df2 <- data.frame(h=h2, l=l2) #V2

h3 <- c("@@", "@@")
l3 <- c("A", "B")
df3 <- data.frame(h=h3, l=l3) #V3

suppRule_Freq <- c(3,0)

outObj <- prepareInput(microDat, filenames=NULL, hierFrames=list(V2=df2,V3=df3,V1=df1), numVar="numVal", suppRule_Freq=suppRule_Freq)
result <- protectTable(outObj, method="HYPERCUBE")
class(result)
summary(result)
}
\keyword{print}