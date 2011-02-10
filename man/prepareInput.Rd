\name{prepareInput}
\alias{prepareInput}
\title{prepareInput}
\description{prepareInput() creates an output object which can be used as \cr
input object in protectTable().}
\usage{prepareInput(dat, filenames=NULL, hierFrames=NULL, freqVar=NULL, numVar=NULL, weightVar=NULL, sampWeightVar=NULL, suppRule_Freq=NULL, suppRule_P=NULL, suppRule_NK=NULL) }
\arguments{
\item{dat}{data.frame or list containing dimensional variables and (optionally) a numeric variable.}
\item{filenames}{if given, vector of filenames of dimensional variables. It is \cr
important that the filename (without extension) equals the variable name which needs\cr
to exist in dat.}
\item{hierFrames}{if given, a list of data-frames given}
\item{freqVar}{if given, the variable name of the variable containing frequency counts for combinations in dat.}
\item{numVar}{if given, the variable name of a numeric variable in dat.}
\item{weightVar}{if given, the variable name of a variable used as weights for HITAS|OPT protecting procedures.}
\item{sampWeightVar}{if given, the variable name of a variable used for weighting purposes.}
\item{suppRule_Freq}{optionally require primary suppression using the freq-rule.}
\item{suppRule_P}{optionally require primary suppression using the p\%-rule.}
\item{suppRule_NK}{optionally require primary suppression using the nk-rule.}
}
\details{
This function generates an output object featuring all kinds of neccessary\cr
information about the dimensional variable under consideration. It should\cr
be noted that the hierarchy-file or the input dataframe given as input\cr
needs to be in correct order. Please have a look at the example-files!}
\value{manipulated data.}

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

outObj <- prepareInput(microDat, filenames=NULL, hierFrames=list(V2=df2,V3=df3,V1=df1), freqVar=NULL, numVar="numVal", weightVar=NULL, sampWeightVar=NULL, suppRule_Freq=suppRule_Freq)
summary(protectTable(outObj, method="HYPERCUBE"))
}
\keyword{methods}