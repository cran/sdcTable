\name{pasteStrVec}
\alias{pasteStrVec}
\title{pasteStrVec}
\description{Pasting string vectors in a specific way.}
\usage{pasteStrVec(strVec, nrVars) }
\arguments{
	\item{strVec}{a character vector.}
	\item{nrVars}{number of columns of the original dataset whose entries should be pasted together}
}
\details{This function is a wrapper function to a c++ function that appends 
vectors. It basically replaces apply(x, 1, paste) calls.}
\value{Manipulated data.}

\author{Bernhard Meindl}
\examples{
### joining 1st and 14, 2nd and 15th element,....
pasteStrVec(letters[1:26], 2)	
}
\keyword{methods}