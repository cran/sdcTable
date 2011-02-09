\name{setBounds}
\alias{setBounds}
\title{setBounds}
\description{sets Bounds needed for the protection of tables with HITAS algorithm.}
\usage{setBounds(outObj, type, v)}
\arguments{
\item{ outObj }{an object derived from calcFullTable().}
\item{ type }{ specifies the type of bounds to be set.\cr Possible values include:\cr'lb' (lower bound known by attacker)\cr'ub' (upper bound known by attacker)\cr'LPL' (lower protection level)\cr'UPL' (upper protection level)\cr'SPL' (sliding protection level)\cr These bounds need to be set for each possible cell.}
\item{ v }{a vector with values for the parameter specified by variable 'type'.}
}
\author{Bernhard Meindl}
\note{the bounds are only used if HITAS is selected as the protection algorithm.\cr
By default, a sliding protection of 1 is set for each cell in calcFullTable()\cr
which means that no primary suppressed cell may be recalculated exactly\cr
after the protection procedure.}
\examples{
	datObj <- get(load(paste(searchpaths()[grep("sdcTable", searchpaths())], "/data/testObj.RData", sep="")))
	# attacker knows that cells are non-negative
	datObj <- setBounds(datObj, type="lb", rep(0, length(datObj$fullTabObj$strID)))
}
\keyword{methods}