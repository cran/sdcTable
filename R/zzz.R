.onLoad <- function(lib, pkg) {
	library.dynam("sdcTable", pkg, lib)
	cat("Package sdcTable has been loaded!\ \n")
}
