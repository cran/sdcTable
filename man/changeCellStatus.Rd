\name{changeCellStatus}
\alias{changeCellStatus}
\title{changeCellStatus}
\description{....}
\usage{changeCellStatus(outObj, varNames, characteristics, rule, codesOrig=TRUE, suppZero=FALSE) }
\arguments{
\item{ outObj }{ an object created by createFullTable() and possible changed by primarySuppression(). }
\item{ varNames }{ a vector of variable names which need to exist in outObj$levelObj. }
\item{ characteristics }{ for the corresponding element in varNames, the characteristic (either original or standardized code) of the cell that needs to be suppressed. }
\item{ rule }{ defines to which status a given cell is changed. Allowed choices are 'u' (mark cell as primary suppressed), 'z' (force publication of cell), 's' (allowed candidate for secondary suppression and 'x' (set cell as secondary suppressed).}
\item{ codesOrig }{ defines if original codes (codesOrig=TRUE) or standardized codes (codesOrig=FALSE) specified in parameter 'characteristics'. }
\item{ suppZero }{ if set to TRUE, a cell is set primary protected even if it has frequency 0. If FALSE, it will not be set protected in this case. }
}
\details{for given characteristics of the dimensional variables spanning the input-table, the index of the cell that needs to be set to primary unsafe is calculated. \cr
The corresponding cell is finally marked as primary unsafe.}
\value{ manipulated data. }
\author{Bernhard Meindl}
\examples{
datObj <- get(load(paste(searchpaths()[grep("sdcTable", searchpaths())], "/data/testObj.RData", sep="")))

# force some cells to be published/not suppressed
varNames <- c("V3","V2","V1","V4")
characteristics <- c("B","m", "010","0011")
rule <- "z"
datObj <- changeCellStatus(datObj, varNames, characteristics, rule, codesOrig=TRUE, suppZero=FALSE)

# set cells to be candidates for secondary suppression
varNames <- c("V1","V2")
characteristics <- c("010","m") # for V3 and V4 the absolute top-level is used
rule <- "s"
datObj <- changeCellStatus(datObj, varNames, characteristics, rule, codesOrig=TRUE, suppZero=FALSE)
}
\keyword{methods}
