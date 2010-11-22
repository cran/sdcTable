## Wrapper function for pasting key-Variables
pasteStrVec <- function(strVec, nrVars) {
	if(length(strVec) %% nrVars != 0)
		stop("Wrong Dimensions!\n")
	else {
		.Call( "myPaste", as.character(strVec), nrVars, PACKAGE = "sdcTable")
	}
}