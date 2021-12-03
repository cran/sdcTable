% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/createJJFormat.R
\name{createJJFormat}
\alias{createJJFormat}
\title{Create input for jj_format}
\usage{
createJJFormat(x)
}
\arguments{
\item{x}{a \linkS4class{sdcProblem} object}
}
\value{
an input suitable for \code{\link[=writeJJFormat]{writeJJFormat()}}
}
\description{
This function transforms a \linkS4class{sdcProblem} object into a list that can
be used as input for \code{\link[=writeJJFormat]{writeJJFormat()}} to write a problem in \code{"JJ-format"} to
disk.
}
\examples{
# setup example problem
# microdata
utils::data("microdata1", package = "sdcTable")

# create hierarchies
dims <- list(
  region = sdcHierarchies::hier_create(root = "Total", nodes = LETTERS[1:4]),
  gender = sdcHierarchies::hier_create(root = "Total", nodes = c("male", "female")))

# create a problem instance
p <- makeProblem(
  data = microdata1,
  dimList = dims,
  numVarInd = "val")

# create suitable input for `writeJJFormat`
inp <- createJJFormat(p); inp

# write files to disk
# frequency table by default
writeJJFormat(
  x = inp,
  path = file.path(tempdir(), "prob_freqs.jj"),
  overwrite = TRUE
)

# or using the numeric variable `val` previously specified
writeJJFormat(
  x = inp,
  tabvar = "val",
  path = file.path(tempdir(), "prob_val.jj"),
  overwrite = TRUE
)
}
\author{
Bernhard Meindl (bernhard.meindl@statistik.gv.at) and
Sapphire Yu Han (y.han@cbs.nl)
}