# update a simple triplet matrix with given col-indices
updateSimpleTripletMatrix <- function(M, i, j, v) {
	if ( !is.simple_triplet_matrix(M) )
		stop("M needs to be a simple triplet matrix!\n")
	
	ind <- which(M$i == i & M$j==j)
	if ( length(ind) > 1 )
		stop("only one element may be updated at one call!\n")
	
	if ( length(ind) == 0 ) {
		if ( v != 0 ) {
			M$i <- append(M$i, i)
			M$j <- append(M$j, j)
			M$v <- append(M$v, v)			
		}
	}
	else {
		if ( v != 0 )
			M$v[ind] <- v		
	}
	M
}

removeFromSimpleTripletMatrix <- function(M, i=NULL, j=NULL) {
	if ( !is.simple_triplet_matrix(M) )
		stop("M needs to be a simple triplet matrix!\n")
	
	if ( !is.null(i) ) {
		indexRows <- which(M$i == i)
		if ( length(indexRows) > 0 ) {
			M$i <- M$i[-indexRows]
			M$j <- M$j[-indexRows]
			M$v <- M$v[-indexRows]
		}		
	}
	if ( !is.null(j) ) {
		indexCols <- which(M$j == j)
		if ( length(indexCols) > 0 ) {
			M$i <- M$i[-indexCols]
			M$j <- M$j[-indexCols]
			M$v <- M$v[-indexCols]
		}		
	}	
	M
}



