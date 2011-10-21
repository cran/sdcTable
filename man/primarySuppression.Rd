\name{primarySuppression}
\alias{primarySuppression}
\title{perform primary suppression in \code{\link{sdcProblem-class}}-objects}
\usage{
  primarySuppression(object, type, ...)
}
\arguments{
  \item{object}{a \code{\link{sdcProblem-class}} object}

  \item{type}{character vector of length 1 defining the
  primary suppression rule. Allowed types are: \itemize{
  \item \code{freq}: apply frequency rule with parameters
  \code{maxN} and \code{allowZeros} \item \code{nk}: apply
  nk-dominance rule with parameters \code{n}, \code{k} and
  \code{numVarInd} \item \code{p}: apply p-percent rule
  with parameters \code{p} and \code{numVarInd} }}

  \item{...}{parameters used in the identification of
  primary sensitive cells. Parameters that can be
  modified|changed are: \itemize{ \item \code{maxN}:
  numeric vector of length 1 used when applying the
  frequency rule. All cells having counts <= \code{maxN}
  are set as primary suppressed. The default value of
  \code{maxN} is 3. \item \code{allowZeros}: logical vector
  of length 1 specifying if empty cells (count==0) should
  be considered sensitive when using the frequency rule.
  The default value of \code{allowZeros} is 'FALSE' so that
  empty cells are not considered primary sensitive by
  default. \item \code{p}: numeric vector of length 1
  specifying parameter \code{p} that is used when applying
  the p-percent rule with default value of 80. \item
  \code{n}: numeric vector of length 1 specifying parameter
  \code{n} that is used when applying the nk-dominance
  rule. Parameter \code{n} is set to 2 by default. \item
  \code{k}: numeric vector of length 1 specifying parameter
  \code{k} that is used when applying the nk-dominance
  rule. Parameter \code{n} is set to 85 by default. \item
  \code{numVarInd}: numeric vector of length 1 specifying
  the index of the numerical variable that should be used
  to identify cells that are dominated by 2 (p-percent
  rule) or n (nk-dominance)-rule. If \code{type} is either
  'nk' or 'p' it is mandatory to specify p. }}
}
\value{
  a \code{\link{sdcProblem-class}} object
}
\description{
  Function \code{\link{primarySuppression}} is used to
  identify and suppress primary sensitive table cells in
  \code{\link{sdcProblem-class}} objects. Argument
  \code{type} allows to select a rule that should be used
  to identify primary sensitive cells. At the moment it is
  possible to identify and suppress sensitive table cells
  using the frequency-rule, the nk-dominance rule and the
  p-percent rule.
}
\note{
  the nk-dominance rule and the p-percent rule can only be
  applied if micro data have been used as input data to
  function \code{\link{makeProblem}}.
}
\examples{
\dontrun{
# load micro data
sp <- searchpaths()
fn <- paste(sp[grep("sdcTable", sp)], "/data/microData1.RData", sep="")
microData <- get(load(fn))

# load problem (as it was created in the example in \\code{\\link{makeProblem}})
fn <- paste(sp[grep("sdcTable", sp)], "/data/problem.RData", sep="")
problem <- get(load(fn))

# we have a look at the table
print(table(microData))

# cell with region=='A' and gender=='female' has 2 units contributing to it
# this cell should be considered senstive!
problem <- primarySuppression(problem, type='freq', maxN=3)

# looking at anonymization states
print(table(getInfo(problem, type='sdcStatus')))

# we see that exactly one cell is primary suppressed (sdcStatus=='u') and
# the remaining cells are possible candidates for secondary suppression ('s')
}
}
\author{
  Bernhard Meindl \email{bernhard.meindl@statistik.gv.at}
}

