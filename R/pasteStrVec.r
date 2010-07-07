## Wrapper function for pasting key-Variables
pasteStrVec <- function(strVec, nrVars) {
	if(length(strVec) %% nrVars != 0)
		stop("Wrong Dimensions!\n")
	else {
		param <- list()
		param$nrKeyVars <- nrVars
		b <- rep("", length=length(strVec)/nrVars)
		out <- .Call( "myPaste", as.character(strVec), b, param, PACKAGE = "sdcTable")$str
	}
	out
}