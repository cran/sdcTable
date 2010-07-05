summary.safeTable <- function(object, ...) {
    nrPrimSupps <- length(which(object$fullData$data$geh=="P"))
    cat("The dataset was protected succesfully using",object$method,"!\n")
    cat("Necessary runs:",object$counter,"\n")
    cat("Total duration:",object$time,"\n")    
    cat("Primary suppressions:", nrPrimSupps , "(", round(100 * (nrPrimSupps /nrow(object$fullData$data)), 2), "%)\n")
    cat("Secondary suppressions:", object$totSupps,"(",round(100*(object$totSupps/nrow(object$fullData$data)),2),"%)\n")
}
