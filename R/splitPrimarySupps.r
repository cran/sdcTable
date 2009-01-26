splitPrimarySupps <- function(fullData) {
    getUpperHier <- function(level, dims) {
        cs <- cumsum(dims)
        if(sum(as.integer(level)) == 0) {
          out <- paste(rep("0", sum(dims)), collapse="")
          return(list(out=out))
        }
        else {
            to <- max(gregexpr("[1-9]", level)[[1]])
            l1 <- min(which(cs >= to))
            from <- cs[l1-1]+1
            len <- to - from + 1
            lev0 <- level;
            substr(lev0, from, to) <- paste(rep("0", length=len), collapse="")
            out <- lev0
        }
        return(out)
    }
    
    fdsub <- fullData$data[fullData$supps2check,]
    ups <- list()
    for (i in 1:fullData$numberindexvars) {
        ups[[i]] <- factor(unlist(lapply(fdsub[,fullData$indexvars[i]], function(x) { getUpperHier(as.character(x), fullData$dims[[i]])  } )) )
    }
    tmp <- as.data.frame(ups)
    colnames(tmp) <- rep("", fullData$numberindexvars)
    str <- "factor(paste(as.character(tmp[,1]), "
    for (i in 2:fullData$numberindexvars) {
        str <- paste(str, "as.character(tmp[,",fullData$indexvars[i],"]), ")
    }
    str <- paste(str, "sep=\"\"))")
    f <- eval(parse(text=str))

    fdsub <- as.data.frame(cbind(fdsub, tmp))
    spl <- split(fdsub, f)
    spl
}
