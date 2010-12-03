.onLoad <- function(lib, pkg) {
	library.dynam("sdcTable", pkg, lib)
	#v <- citation("sdcTable")$note
	#version <- substr(v, nchar(v)-4, nchar(v))
	#cat("Package sdcTable", version, "has been loaded!\n")
	cat("Package sdcTable 0.0.15 has been loaded!\n")
}
