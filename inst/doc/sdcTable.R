## ---- include=FALSE-----------------------------------------------------------
library(sdcTable)

## ---- eval=TRUE---------------------------------------------------------------
library(sdcTable)
packageVersion("sdcTable")

## ---- include=FALSE-----------------------------------------------------------
# daten laden
microData <- readRDS("microdat.rds")
aggregatedData <- readRDS("aggdat.rds")
completeData <- readRDS("completedat.rds")

h1 <- sdcHierarchies::hier_create("Tot", nodes = LETTERS[1:4])
h1 <- sdcHierarchies::hier_add(h1, "B", nodes = paste0("B", letters[1:3]))
h2 <- sdcHierarchies::hier_create("Tot", nodes = c("m", "w"))
h3 <- sdcHierarchies::hier_create("Tot", nodes = letters[1:6])
dimList <- list(
  V1 = sdcHierarchies::hier_convert(h1, as = "df"),
  V2 = sdcHierarchies::hier_convert(h2, as = "df"),
  V3 = sdcHierarchies::hier_convert(h3, as = "df")
)

#resOPT <- workDat$resOPT
#completeData <- workDat$resOPT@finalData
#completeData <- completeData[, -ncol(completeData)]

## ---- echo=FALSE--------------------------------------------------------------
print(head(microData), row.names = FALSE)

## ---- echo=FALSE--------------------------------------------------------------
lev.V1 <- as.character(sort(unique(microData$V1)))
lev.V2 <- as.character(sort(unique(microData$V2)))
lev.V3 <- as.character(sort(unique(microData$V3)))

## ---- echo=FALSE, comment=""--------------------------------------------------
cat(paste(shQuote(lev.V1), collapse = ", "))

## ---- echo=FALSE, comment=""--------------------------------------------------
cat(paste(shQuote(lev.V2), collapse = ", "))

## ---- echo=FALSE, comment=""--------------------------------------------------
cat(paste(shQuote(lev.V3), collapse = ", "))

## -----------------------------------------------------------------------------
print(tail(completeData))

## ---- echo=FALSE--------------------------------------------------------------
levComp.V1 <- dimList$V1$name
levComp.V2 <- dimList$V2$name
levComp.V3 <- dimList$V3$name

## ---- echo=FALSE, comment=""--------------------------------------------------
cat(paste(shQuote(levComp.V1), collapse = ", "))

## ---- echo=FALSE, comment=""--------------------------------------------------
cat(paste(shQuote(levComp.V2), collapse = ", "))

## ---- echo=FALSE, comment=""--------------------------------------------------
cat(paste(shQuote(levComp.V3), collapse = ", "))

## ---- echo=FALSE--------------------------------------------------------------
x <- completeData[nrow(completeData), ]

## ---- echo=FALSE, comment=""--------------------------------------------------
subTots.V1 <- setdiff(levComp.V1, lev.V1)
cat(paste(shQuote(subTots.V1), collapse = ", "))

## ---- echo=FALSE, comment=""--------------------------------------------------
subTots.V2 <- setdiff(levComp.V2, lev.V2)
cat(paste(shQuote(subTots.V2), collapse = ", "))

## ---- echo=FALSE, comment=""--------------------------------------------------
subTots.V3 <- setdiff(levComp.V3, lev.V3)
cat(paste(shQuote(subTots.V3), collapse = ", "))


## -----------------------------------------------------------------------------
dimV1 <- matrix(nrow = 0, ncol = 2)
dimV1 <- rbind(dimV1, c("@", "Tot"))
print(dimV1)

## -----------------------------------------------------------------------------
mat <- matrix(nrow = 4, ncol = 2)
mat[, 1] <- rep("@@", 4)
mat[, 2] <- LETTERS[1:4]
dimV1 <- rbind(dimV1, mat)
print(dimV1)

## -----------------------------------------------------------------------------
mat <- matrix(nrow = 3, ncol = 2)
mat[, 1] <- rep("@@@", 3)
mat[, 2] <- c("Ba", "Bb", "Bc")

dimV1 <- rbind(dimV1, mat)
print(dimV1)

## -----------------------------------------------------------------------------
dimV1 <- dimV1[c(1:3,6:8, 4:5),]
print(dimV1, row.names = FALSE)

## -----------------------------------------------------------------------------
dimV1 <- sdcHierarchies::hier_create(root = "Tot", nodes = LETTERS[1:4])
dimV1 <- sdcHierarchies::hier_add(dimV1, root = "B", nodes = c("Ba","Bb","Bc"))
sdcHierarchies::hier_display(dimV1)

## -----------------------------------------------------------------------------
dimV2 <- sdcHierarchies::hier_create(root = "Tot", nodes = c("m", "w"))
sdcHierarchies::hier_display(dimV2)

## -----------------------------------------------------------------------------
dimV3 <- sdcHierarchies::hier_create(root = "Tot", nodes = letters[1:6])
sdcHierarchies::hier_display(dimV3)

## -----------------------------------------------------------------------------
dimList <- list(V1 = dimV1, V2 = dimV2, V3 = dimV3)
prob.microDat <- makeProblem(
	data = microData,
	dimList = dimList,
	dimVarInd = 1:3,
	freqVarInd = NULL,
	numVarInd = 4:5,
	weightInd = NULL,
	sampWeightInd = NULL)

## -----------------------------------------------------------------------------
### problem from complete data ###
dimList <- list(V1 = dimV1, V2 = dimV2, V3 = dimV3)
prob.completeDat <- makeProblem(
	data = completeData,
	dimList = dimList,
	dimVarInd = 1:3,
	freqVarInd = 4,
	numVarInd = 5:6,
	weightInd = NULL,
	sampWeightInd = NULL)

## -----------------------------------------------------------------------------
all(c(class(prob.microDat), class(prob.completeDat)) == "sdcProblem")

## -----------------------------------------------------------------------------
counts1 <- getInfo(prob.completeDat, type = "freq")
counts2 <- getInfo(prob.microDat, type = "freq")
all(counts1 == counts2)

## -----------------------------------------------------------------------------
prob.completeDat <- primarySuppression(prob.completeDat, type = "freq", maxN = 10)

## -----------------------------------------------------------------------------
print(table(getInfo(prob.completeDat, type = "sdcStatus")))
summary(prob.completeDat)

## ---- echo=FALSE--------------------------------------------------------------
nrPrimSupps <- length(which(getInfo(prob.completeDat, type = "sdcStatus") == "u"))

## ---- cache=TRUE--------------------------------------------------------------
resHITAS <- protectTable(prob.completeDat, method = "HITAS")
resOPT <- protectTable(prob.completeDat, method = "OPT")
resHYPER <- protectTable(prob.completeDat, method = "HYPERCUBE")
resSIMPLE <- protectTable(prob.completeDat, method = "SIMPLEHEURISTIC")

## -----------------------------------------------------------------------------
finalData <- getInfo(resOPT, type = "finalData")
print(head(finalData))

## -----------------------------------------------------------------------------
summary(resOPT)

