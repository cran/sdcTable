\docType{methods}
\name{set.cutList}
\alias{set.cutList}
\alias{set.cutList,cutList,character,list-method}
\title{modify \code{cutList}-objects depending on argument \code{type}}
\arguments{
  \item{object}{an object of class \code{cutList}}

  \item{type}{a character vector of length 1 defining what
  to calculate|return|modify. Allowed types are:} \itemize{
  \item addCompleteConstraint: add a constraint to argument
  \code{object} \item removeCompleteConstraint: remove a
  constraint from argument \code{object}}

  \item{input}{a list depending on argument \code{type}.}
  \itemize{ \item type==addCompleteConstraint: a list of
  length 1 \itemize{ \item first element: an object of
  class \code{cutList} with exactly one constraint } \item
  type==removeCompleteConstraint: a list of length 1
  \itemize{ \item first element: numeric vector of length 1
  specifying the index of the constraint that should be
  removed }}
}
\value{
  an object of class \code{cutList}
}
\description{
  modify \code{cutList}-objects depending on argument
  \code{type}
}
\note{
  internal function
}
\author{
  Bernhard Meindl \email{bernhard.meindl@statistik.gv.at}
}

