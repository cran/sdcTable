\docType{methods}
\name{set.dataObj}
\alias{set.dataObj}
\alias{set.dataObj,dataObj,character,listOrNULL-method}
\title{modify \code{dataObj}-objects depending on argument \code{type}}
\arguments{
  \item{object}{an object of class \code{dataObj}}

  \item{type}{a character vector of length 1 defining what
  to calculate|return|modify. Allowed types are:} \itemize{
  \item rawData: set slot 'rawData' of argument
  \code{object}}

  \item{input}{a list depending on argument \code{type}.}
  \itemize{ \item type==rawData: a list containing raw
  data}
}
\value{
  an object of class \code{dataObj}
}
\description{
  modify \code{dataObj}-objects depending on argument
  \code{type}
}
\note{
  internal function
}
\author{
  Bernhard Meindl \email{bernhard.meindl@statistik.gv.at}
}

