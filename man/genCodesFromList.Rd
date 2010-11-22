\name{genCodesFromList}
\alias{genCodesFromList}
\title{genCodesFromList}
\description{Generate standardized codes needed from a hierarchy-file or a\cr
input dataframe.}
\usage{genCodesFromList(file=NULL, dataframe=NULL, hierIndex="@", varname, minimal=FALSE) }
\arguments{
\item{file}{path to hierarchy-file.}
\item{dataframe}{dataframe with 2 columns representing the hierarchy levels and the characteristics.}
\item{hierIndex}{character used to specify different levels in hierarchy-file}
\item{varname}{variable-name to be standardized}
\item{minimal}{if TRUE, all (sub)totals will be removed}
}
\details{
This function generates an output objects similar to the 
input-hierarchy-file with an additional column specifiying 
the standardized codes. It should be noted that the hierarchy-file
or the input dataframe given as input need to be in correct order. 
Please have a look at the example-files!}
\value{Manipulated data.}

\author{Bernhard Meindl}
\examples{
### calculate standardized code from example hierarchy-file
file <-  paste(searchpaths()[grep("sdcTable", searchpaths())], "/etc/exampleHier.hcr", sep="")
dataF <- read.table(file, sep=";", dec=".", colClasses="character")	

# directly from file
all.codes <- genCodesFromList(file=file, hierIndex="@", varname="someVar")

# from the table
all.codes2 <- genCodesFromList(dataframe=dataF, hierIndex="@", varname="someVar")

minimal.codes <- genCodesFromList(file=file, hierIndex="@", varname="someVar", minimal=TRUE)
print(str(all.codes))
codesOrig <- minimal.codes$codesOrig
codesStandard <- minimal.codes$codesStandard
levelStructure <- minimal.codes$levelStructure
print(codesOrig)
print(codesStandard)
print(levelStructure)
}
\keyword{methods}