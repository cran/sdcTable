\name{createFullData}
\alias{createFullData}
\title{createFullData}
\description{Function to generate an object of class fullData which can then be used as input object for several other functions in sdcTable such as protectTable.}
\usage{createFullData (minimalData, indexvars, l, suppVals=FALSE, suppLimit=NULL, suppZeros=NULL) }
\arguments{
  \item{minimalData}{ a data.frame containing a column for each hierarchical variable in a specific coding and a column for the corresponding numerical values. The values need to be integer at the moment.}
  \item{indexvars}{a vector containing the position of the hierarchical variables within minimalData }
  \item{l}{a list which element i defining the hierarchical structure of the i-th hierarchical variable. }
  \item{suppVals}{should during the generation of the complete data set values be suppressed? }
  \item{suppLimit}{if suppVals is TRUE, what is the threshold for primary suppressing values? () }
  \item{suppZeros}{binary: should values of empty cells be suppressed as well? }
}
\details{Have a look at the link given below.}
\value{object of class fullData}
\author{Bernhard Meindl}
\note{fix me: LOTS! more primary suppression rules should be implemented; what to do with real numbers? simplify the generation of object of class fullData? }
\examples{
	\dontrun{
		# create minimalData (simple 2-dimensional Example)
		# first hierarchical variable: sex: we have {total}-> {male, women}
		sex <- c("00", "01", "02")
		
		# second hierarchical variable: a total of 13 age-groups including the total: we have {total}<- {group1-group12} 
		age <- c("000", "001", "002", "003", "004", "005", "006", "007", "008", "009", "010", "011", "012")
		
		# for the minimal dataset we don't need the totals: we generate each combination of sublevels of the hierarchical variables
		minimalData <- expand.grid(sex[-1], age[-1])
		
		# the column with values is added (just random numbers in this case)
		minimalData$value <- rpois(nrow(minimalData), 7)
		
		# position of the hierarchical variables in minimalData
		indexvars <- c(1,2)		
		
		# level-structure of the hierarchical variables:
		# the first level (sex) consists of a total of 2 levels (this is total and male/women). Since there are only 2 characteristics in the sub-level, we need to add 1.
		# the first level (age)consists of a total of 12 levels (this is the total and 12 age groups). Since we need 2 digits for the 12 age groups (01, 02,...,12)for this level, the second component of list[[2]] needs to be 2.
		l <- list()
		l[[1]] <- c(1, 1)
		l[[2]] <- c(1, 2)
		
		# we generate an object of class fullData which we can then use as input object for other functions and primary suppress all values less equal than 4, but don't mark emty cells as primary suppressed.
		fullDat <- createFullData (minimalData, indexvars, l, suppVals=TRUE, suppLimit=4, suppZeros=FALSE)
		class(fullDat)
		
		#result <- protectTable(fullDat, method="GHMITER")
		summary(result)
		
		# result$fullData$data is the protected dataset containing all totals. The user needs then to provide a function to recode the hierarchical variables.
	}
}
\keyword{methods}