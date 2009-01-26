## calcSubset is a function to draw a specific subtable from the entire data-structure
## fullData: the entire data structure as given in Quaderverfahren() and qf()
## group: list-element from output of splitPrimarySupps()
calcSubset <- function(fullData, group) {
    # getLevels() and getSameHier() are only used as auxiliary functions
    getLevels <- function(sub, indexvars, j) {
        tmp <- sub[j, indexvars]
        out <- NULL
        for (i in 1:length(indexvars)) {
            out <- c(out, as.character(sub[j,][[indexvars[i]]]))
        }
        return(out)
    }
    getSameHier <- function(level, dim) {
        out <- list()
        cs <- cumsum(dim)
        if(sum(as.integer(level)) == 0) {
          out[[1]] <- -1
          out[[2]] <- (cs[2]+1):cs[length(cs)]
        }
        else {
            to <- max(gregexpr("[1-9]", level)[[1]])
            l1 <- min(which(cs >= to))

            out[[1]] <- 1:cs[(l1-1)]
            if(l1 < length(cs)) {
                out[[2]] <- ((cs[l1])+1):cs[length(cs)]
            }
            else {
                out[[2]] <- -1
            }
        }
        return(out)
    }

    subtab <- fullData
    index <- getLevels(group, fullData$indexvars, 1)
    for (i in 1:fullData$numberindexvars) {
        s <- getSameHier(index[i], fullData$dims[[i]])
        # data subset in which characters from:to are equal to the current level
        if(s[[1]][1] != -1) {
            from <- s[[1]][1]
            to <-  s[[1]][length(s[[1]])]
            subtab$data <- subset(subtab$data, substr(subtab$data[,fullData$indexvars[i]], from, to) == substr(index[i], from, to))
        }
        # sometimes we need to check for "00" at the end
        if(s[[2]][1] != -1 ) {
            from <- s[[2]][1]
            to <-  s[[2]][length(s[[2]])]
            subtab$data <- subset(subtab$data, substr(subtab$data[,fullData$indexvars[i]], from, to) == paste(rep("0", to-from+1), collapse=""))
        }
    }
    return(subtab)
}