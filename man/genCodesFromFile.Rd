\name{genCodesFromFile}
\alias{genCodesFromFile}
\title{genCodesFromFile}
\description{Generate standardized codes needed from a hierarchy-file.}
\usage{genCodesFromFile(file, hierIndex="@", varname, minimal=FALSE) }
\arguments{
  \item{file}{path to hierarchy-file.}
  \item{hierIndex}{character used to specify different levels in hierarchy-file}
  \item{varname}{variable-name to be standardized}
  \item{minimal}{if TRUE, all (sub)totals will be removed}
}
\details{This function generates an output objects similar to the input-hierarchy-file with an additional column specifiying the standardized codes}
\value{Manipulated data.}

\author{Bernhard Meindl}
\examples{
	\dontrun{
		### calculate standardized code from example hierarchy-file
		file <-  paste(searchpaths()[grep("sdcTable", searchpaths())], "/etc/exampleHier.hcr", sep="")
		all.codes <- genCodesFromFile(file=file, hierIndex="@", varname="someVar")
		minimal.codes <- genCodesFromFile(file=file, hierIndex="@", varname="someVar", minimal=TRUE)
		print(str(all.codes))
		codesOrig <- minimal.codes$codesOrig
		codesStandard <- minimal.codes$codesStandard
		levelStructure <- minimal.codes$levelStructure
		print(codesOrig)
		print(codesStandard)
		print(levelStructure)
	}	
}
\keyword{methods}