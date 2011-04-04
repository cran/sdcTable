\name{calcDimInfos}
\alias{calcDimInfos}
\title{calcDimInfos}
\description{calcDimInfos() calculates all necessary information about a \cr
	dimensional variable, eg. this function calculates standardized codes\cr
	from an input-file or an input-object as well as the position of the\cr
	dimensional variable within the data-set that needs to be protected,\cr
	the level-structure or a complete listing of all (sub)-levels.}
\usage{calcDimInfos(inputDat, file=NULL, dataframe=NULL, vName) }
\arguments{
\item{inputDat}{an input-data object (microdata or data-frame)}
\item{file}{path to a dimension (hierarchy)-file.}
\item{dataframe}{dataframe with 2 columns representing the levels and \cr 
the characteristics or a dimension.}
\item{vName}{variable-name of the dimension that is dealt with in inputDat}
}
\details{
This function generates an output object featuring all kinds of neccessary\cr
information about the dimensional variable under consideration. It should\cr
be noted that the hierarchy-file or the input dataframe given as input\cr
needs to be in correct order. Please have a look at the example-files!}
\value{manipulated data.}

\author{Bernhard Meindl}
\examples{
n <- 15
dat <- data.frame(
	sex=sample(c("Male","Female"), n, replace=TRUE),
	age=sample(paste("Age Group",1:12), n, replace=TRUE)
)

### calculate standardized code from example hierarchy-file
# directly from a hierarchy-file
file.Sex <-  paste(searchpaths()[grep("sdcTable", searchpaths())], "/etc/exampleSex.hcr", sep="")
dim.Sex <- calcDimInfos(dat, file=file.Sex, vName="sex")
print(dim.Sex)

# from a data-frame
dataAge <- read.table(paste(searchpaths()[grep("sdcTable", searchpaths())], "/etc/exampleAge.hcr", sep=""), sep=";", dec=".", colClasses="character")	
dim.Age <- calcDimInfos(dat,dataframe=dataAge, vName="age")
print(dim.Age)

# duplicate levels ("BroadAgeGroup 1" is identical to Total Age"
# thus it is listed as list-element "dups" in dim.Age2
dataAge$V1[2:nrow(dataAge)] <- "@@@" 
dataAge <- rbind(dataAge[1,], c("@@", "broad Age Group 1"), dataAge[2:nrow(dataAge),])
dim.Age2 <- calcDimInfos(dat,dataframe=dataAge, vName="age")
print(dim.Age2)
}
\keyword{methods}