.onLoad <- function(lib, pkg) {
	library.dynam("sdcTable", pkg, lib)
	
	myGlobalEnv <- indCplex <- indSymphony <- indLpSolveAPI <- indSnowfall <- NULL
	
	#suppressWarnings(indCplex <- require(Rcplex, quietly=TRUE))
	suppressWarnings(indSymphony <- require(Rsymphony, quietly=TRUE))
	suppressWarnings(indLpSolveAPI <- require(lpSolveAPI, quietly=TRUE))
	suppressWarnings(indSnowfall <- require(snowfall, quietly=TRUE))
	
	myGlobalEnv <- new.env()
	attach(myGlobalEnv)	
	indCplex <- FALSE # temp: remove Cplex support
	assign("indCplex", indCplex, pos=which(search()=="myGlobalEnv"))
	assign("indSymphony", indSymphony, pos=which(search()=="myGlobalEnv"))
	assign("indLpSolveAPI", indLpSolveAPI, pos=which(search()=="myGlobalEnv"))
	assign("indSnowfall", indSnowfall, pos=which(search()=="myGlobalEnv"))
	
	#cat("Is package 'Rcplex' available? --> ")
	#if ( indCplex==TRUE )
	#	cat("Yes (package loaded)")
	#else cat("No")
	#cat("\n")
	
	cat("Is package 'RSymphony' available? --> ")
	if ( indSymphony==TRUE )
		cat("Yes (package loaded)")
	else cat("No")	
	cat("\n")
	
	cat("Is package 'lpSolveAPI' available? --> ")
	if ( indLpSolveAPI==TRUE )
		cat("Yes (package loaded)")
	else cat("No")	
	cat("\n")
	
	cat("Is package 'snowfall' available? --> ")
	if ( indSnowfall==TRUE )
		cat("Yes (package loaded)")
	else cat("No")	
	cat("\n\n")		
		
	cat("Package sdcTable 0.6.4 has been loaded!\n")
	cat("Note: RCplex support not yet available. Please manually install and load this package if you need to.\n")
}
