\docType{methods}
\name{set.dimInfo}
\alias{set.dimInfo}
\alias{set.dimInfo,dimInfo,character,character-method}
\title{modify \code{dimInfo}-objects depending on argument \code{type}}
\arguments{
  \item{object}{an object of class \code{dimInfo}}

  \item{type}{a character vector of length 1 defining what
  to calculate|return|modify. Allowed types are:} \itemize{
  \item strID: set slot 'strID' of argument \code{object}}

  \item{input}{a list depending on argument \code{type}.}
  \itemize{ \item type==strID: a character vector
  containing ID's}
}
\value{
  an object of class \code{dimInfo}
}
\description{
  modify \code{dimInfo}-objects depending on argument
  \code{type}
}
\note{
  internal function
}
\author{
  Bernhard Meindl \email{bernhard.meindl@statistik.gv.at}
}

