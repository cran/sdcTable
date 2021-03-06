% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics_problemInstance.r,
%   R/methods_class_problemInstance.r
\docType{methods}
\name{calc.problemInstance}
\alias{calc.problemInstance}
\alias{calc.problemInstance,problemInstance,character,list-method}
\title{perform calculations on \code{problemInstance}-objects depending on argument \code{type}}
\usage{
calc.problemInstance(object, type, input)

\S4method{calc.problemInstance}{problemInstance,character,list}(object, type, input)
}
\arguments{
\item{object}{an object of class \code{problemInstance}}

\item{type}{a character vector of length 1 defining what to calculate|return|modify. Allowed types are:}
\itemize{
\item makeMasterProblem: create the master problem that is the core of the secondary cell suppression problem
\item isProtectedSolution: check if a solution violates any required (upper|lower|sliding) protection levels}

\item{input}{a list depending on argument \code{type}.}

\itemize{
\item type==makeMasterProblem: input is not used (empty list)
\item type==isProtectedSolution: input is a list of length 2 with elements 'input1' and 'input2'
\itemize{
\item element 'input1': numeric vector of calculated known lower cell bounds (from attacker's problem)
\item element 'input2': numeric vector of known upper cell bounds (from attacker's problem) }}
}
\value{
information from objects of class \code{problemInstance} depending on argument \code{type}
\itemize{
\item an object of class \code{linProb} if argument \code{type} matches 'makeMasterProblem'
\item logical vector of length 1 if argument \code{type} matches 'isProtectedSolution' with TRUE if all primary suppressed cells are adequately protected, FALSE otherwise }
}
\description{
perform calculations on \code{problemInstance}-objects depending on argument \code{type}
}
\note{
internal function
}
\author{
Bernhard Meindl \email{bernhard.meindl@statistik.gv.at}
}
\keyword{internal}
