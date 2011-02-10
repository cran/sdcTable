summary.safeTable <- function(object, ...) {
    nrPrimSupps <- sum(object$nrPrimarySupps)
	nrSecondarySupps <- sum(object$nrSecondSupps)
	if ( object$method == "HITAS" )
    	cat("The dataset was protected succesfully using 'HITAS' ( solver =",object$solver,")!\n")
	if ( object$method == "HYPERCUBE" )
		cat("The dataset was protected succesfully using 'HYPERCUBE'!\n")
	#cat("Necessary runs:",object$counter,"\n")
    cat("Total duration:",object$time,"\n")    
    cat("Primary suppressions:", nrPrimSupps , "(", round(100 * (nrPrimSupps /length(object$outObj$strID)), 2), "%)\n")
    cat("Secondary suppressions:", nrSecondarySupps,"(",round(100*(nrSecondarySupps/length(object$outObj$strID)),2),"%)\n")
	cat("Total suppressions:", nrPrimSupps+nrSecondarySupps,"(",round(100*((nrPrimSupps+nrSecondarySupps)/length(object$outObj$strID)),2),"%)\n")

	
	#setTxtProgressBar(a, round(100 * (nrPrimSupps/length(object$outObj$strID)), 2))
	#cat("bla1\n"); close(a)
	#setTxtProgressBar(b, round(100 * (nrSecondarySupps/length(object$outObj$strID)), 2))
	#cat("bla2\n"); close(b)
	#setTxtProgressBar(c,  round(100 * ((nrPrimSupps + nrSecondarySupps)/length(object$outObj$strID)), 2))
	#cat("bla3\n"); close(c)
	
	#va <- round(100 * (nrPrimSupps/length(object$outObj$strID)), 2)
	#vb <- round(100 * (nrSecondarySupps/length(object$outObj$strID)), 2)
	#vc <- va+vb
	#a <- txtProgressBar(0, 100, va, style=3); cat(" primary suppressions\n")
	#b <- txtProgressBar(0, 100, vb, style=3); cat(" secondary suppressions\n")
	#c <- txtProgressBar(0, 100, vc, style=3); cat(" total suppressions\n")
}
