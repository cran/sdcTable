setBounds <- function(outObj, type, v) {
	if ( !type %in% c("lb", "ub", "LPL", "UPL", "SPL") )
		stop("please check input-parameter 'type'!\n")
	
	if ( length(v) != length(outObj$fullTabObj$strID) )
		stop("please check the length of input 'v'!\n")
	
	if( type == "lb")
		outObj$fullTabObj$lb <- v
	if( type == "ub")
		outObj$fullTabObj$ub <- v	
	if( type == "LPL" ) 
		outObj$fullTabObj$LPL <- v	
	if( type == "UPL" ) 
		outObj$fullTabObj$UPL <- v	
	if( type == "SPL" ) 
		outObj$fullTabObj$SPL <- v		
	outObj
}