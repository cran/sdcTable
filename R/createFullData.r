createFullData <- function(minimalData, indexvars, l, suppVals=FALSE, suppLimit=NULL, suppZeros=NULL) {
    genCompleteHier <- function(sub, hierUniques, l, indexvar) {
        # auxiliary function that calculates the entire hierarchical structure of a given parameter value
        completeHier <- function(str, dim) {
          tmp <- str
          hier <- str
          for (i in length(dim):1) {
            if(dim[i] == 1) {
              from <- sum(dim[1:i])
              to <- from
            }
            else {
              to <- sum(dim[1:i])
              from <- to - dim[i] + 1
            }

            substr(tmp, from, to) <- paste(rep("0", (to-from)+1), collapse="")
            hier <- c(hier, tmp)
          }
          hier <- c(paste(rep("0", sum(dim)), collapse=""), hier)
          return(sort(unique(hier)))
        }
        # sub: raw data without hierachies (lowest levels only)
        # hierUnique: vector of all uniques values of a hierarchical variable
        # l <- how is the hierarchie built up?
        # indexvar: columns of index-variables
        reduceHiers <- function(sub, hierUnique, l, indexvar) {
            sub[,indexvar] <- as.character(sub[,indexvar])
            cs <- cumsum(l)
            index <- length(cs)-1           # set index at penultimate level
            dups <- NULL
            out <- list()
            z <- 1
            Indikator <- TRUE
            while (Indikator!=FALSE) {
                spl <- split(hierUnique, substr(hierUnique, 1, cs[index]))
                dups <- NULL
                for (i in 1:length(spl)) {
                    # the second part of the if-clause means that e.g 0310 must not be added if 0311 does not exist!
                    if(length(spl[[i]]) == 1 & sum( as.integer(substr(as.character(spl[[i]][1]), cs[index]+1, sum(l))) != 0)) {
                        dups <- c(dups, spl[[i]])
                    }
                }

                ifelse(length(dups) == 0, Indikator <- FALSE, Indikator <- TRUE)
                if (Indikator == TRUE) {
                    ind <- which(hierUnique %in% dups)
                    hierUnique <- hierUnique[-ind]
                    out[[z]] <- dups  # save: these characteristics need to be added later
                    z <- z + 1

                    from <- cs[index]+1
                    to <- from + l[index+1] -1

                    substr(dups,from,to) <- paste(rep("0", length=(to-from+1)), collapse="")
                    hierUnique <- c(hierUnique, dups)
                    # adjust sub
                    ind <- which(substr(sub[,indexvar],1,cs[index]) %in% substr(dups,1,cs[index]))
                    substr(sub[,indexvar][ind],from,to) <- paste(rep("0", length=(to-from+1)), collapse="")
                    index <- index - 1
                }
            }
            return(list(sub=sub, hierUnique=sort(hierUnique), toAdd=out))
        }

        # 1) remove redundant hierarchies using reduceHiers()
        reduceUniques <- reduceHiers(sub, hierUniques, l, indexvar)
        sub <- reduceUniques$sub
        hierUniques <- reduceUniques$hierUnique
        hierToAdd <- reduceUniques$toAdd

        # 2) create complete structure of hierarchy using completeHier()
        hier <- NULL
        for (i in 1:length(hierUniques)) {
            hier <- c(hier, completeHier(hierUniques[i], l))
        }
        hier <- sort(unique(hier))

        # return results
        return(list(completeHier=hier, sub=sub, hierToAdd=hierToAdd))
    }

    ## Create the complete hierarchie for a given subset of data)
    ## Input-parameters:
    # sub:          the data.frame containing the most detailed structure of the hierarchies (
    #               (from this data.set, the higher hierarchies can be computed)
    # dims:         list of unique characteristics for each of the hierarchical variables
    # l             list of indices of the hierarchical variables. tells, how many characters are needed for each level
    # indexvars     a vector of index-variables
    createFullHier <- function(sub, dims, l, indexvars) {
        ## Sub-Function for create.full.hier()
        # Based on the necessary subset of data, the full data.set is filled with values from sub
        ## Input-Parameters:
        # sub:        the smallest possible subset of data
        # structure:  the full dataset to be filled with values from sub
        # indexvars:  the columns of the index-variables
        fillStructure <- function(sub, structure, indexvars) {
            #  we fill the complete dataset with values from "sub"
            for (i in 1:nrow(sub)) {
                x <- as.matrix(sub[i,indexvars])[indexvars]

                str <- paste("which(structure[,", indexvars[1], "] ==  \"", x[1], "\"", sep="")
                if(length(x) >= 2) {
                    for (j in 2:length(x)) {
                        str <- paste(str, " & structure[,", indexvars[j], "] == \"", as.character(x[j]), "\"", sep="")
                    }
                }
                str <- paste(str, ")", sep="")
                element <- eval(parse(text=str))
                structure[element,"val"] <- sub[i, "val"]
            }
            return(structure)
        }

        # fill.vals() is a function which calculates the values needed for a given factor based on it's subfactors(for a given index.variable and subset of data)
        ## Input-Parameters:
        # str:          the factor under consideration
        # sub:          a dataset containing the values which need to be summed up to be correct for "factor")
        # indexvars:    a vector of indices belonging to the hierarchical variables
        # akt.ind:      the hierarchical variable under consideration
        # structure:    the data.set that needs to be filled with (summed) values from sub
        fillVals <- function(str, l, sub, indexvars, akt.ind, structure) {
                ## Subfunction that is called from fill.vals()
                # Returns the first index from factor "str" (with length l) that is != 0
            lastNotNull <- function(str, l) {
                if(as.integer(str) != 0) {
                    ind <- NULL
                    level <- 1
                    u <- 0
                    while (is.null(ind)) {
                        u <- u + l[level]
                        diff <- l[level] - 1
                        sub  <- as.numeric(substr(str, u-diff, u))

                        if( (level+1) <= length(l) ) {
                            u2 <- u+ l[level+1]
                            diff2 <- l[level+1]-1
                            sub2 <- as.numeric(substr(str, u2-diff2, u2))
                            ifelse (sub > 0 & sub2==0, ind <- 1, level <- level +1 )
                            tmp <- substr(str, 1, sum(l[1:level]))
                        }
                        else {
                            ifelse (sub > 0, ind <- 1, level <- level +1 )
                            tmp <- substr(str, 1, sum(l[1:level]))
                        }
                    }
                }
                # special case: top-hierarchy (-> we report negative values and make a switch in function get.subset())
                else {
                    sub <- -1
                    level <- -1
                    tmp <- paste(rep("0", sum(l)), collapse="")
                }
                return(list(level=level, value=sub, l=l, tmp=tmp))
            }

            ## subfunction that is called from fillVals()
            # Returns a subset of data given an index for the hierarchical variable and an object of the function "firstNotNull()
            getSubset <- function(x, sub, index) {
                s <- subset(sub, sub$val != -1)
                level <- x$level
                subval <- x$value
                tmp <- x$tmp
                l <- x$l
                ind <- index
                # top-level
                if (level==-1 & subval == -1) {
                    ind2 <- sum(l[1:2])
                    if(length(l) > 2) {
                        sub.new <- subset(s, substr(as.character(s[,ind]), l[1]+1, ind2) !="0" & as.integer(substr(s[,ind], sum(l[1:3]), sum(l))) == 0)
                    }
                    else {
                        sub.new <- subset(s, substr(as.character(s[,ind]), l[1]+1, ind2) !="0")
                    }
                }
                # not the top-level!
                else {
                    sub.new <- subset(s, substr(as.character(s[,ind]), 1, nchar(tmp)) == tmp)

                    # if there are subsublevels <- we have to exclude them....
                    if((level+2) <= length(l)) {
                        unten <- cumsum(l)[level+1]+1
                        oben <- unten + (l[level+2]-1)
                        sub.new <- subset(sub.new, as.integer(substr(as.character(sub.new[,ind]), unten, oben)) == 0)
                    }
                }
                return(sub.new)
            }

            ## Subfunction that is called from fillVals()
            # splits the data.set returned by get.subset() according to the factors of the other indexvars
            splitSubs <- function(new.subs, indexvars, akt.var) {
                index <- indexvars[-akt.var]
                str <- paste("paste(as.character(new.subs[,", index[1], "])" , sep="")
                l <- length(index)

                if (l > 1) {
                    for (i in 2:l) {
                        str <- paste(str, ", as.character(new.subs[,", index[i], "])" , sep="")
                    }
                }
                str <- paste(str, ", sep='')", sep="")
                the_call <- paste("factor(", str, ")", sep="")
                new.subs$factor <- eval(parse(text=the_call))

                split <- split(new.subs$val, f=new.subs$factor)
                return(list(split=split, sub=new.subs, indexvars=indexvars, akt.var=akt.var))
            }

            ## Subfunction that is called from fillVals()
            # takes an object from split.subs() as input and adds corresponding values to the complete dataset (structure)
            sumSplitSubs <- function(split, structure, string) {
                indexvars <- split$indexvars
                akt.var <- split$akt.var
                index <- indexvars[-akt.var]

                new.subs <- split$sub
                split.new <- split$split

                for (j in 1:length(split.new)) {
                    fac <- names(split.new[j])
                    sum <- sum(split.new[[j]])

                    new.subset <- subset(new.subs, new.subs$factor==fac)[1,indexvars[-akt.var]]

                    # we calculate the indices we can insert in structure
                    str <- paste("row.names(subset(structure, structure[,", akt.var, "] == '", as.character(string), "'",  sep="")
                    for (i in 2:length(indexvars)) {
                        if(is.null(nrow(new.subset))) {
                            str <- paste(str, " & structure[, ", index[i-1], "] == as.character('", new.subset, "')", sep="")
                        }
                        else {
                            str <- paste(str, " & structure[, ", index[i-1], "] == as.character('", new.subset[1,i-1], "')", sep="")
                        }
                    }
                    the_call <- paste(str, "))", sep="")
                    index.structure <- as.integer(eval(parse(text=the_call)))

                    # we add the values to "structure" (with the indices calculated above)
                    structure$val[index.structure] <- sum
                }
                return(structure)
            }

            ind <- lastNotNull(str, l)
            newSubs <- getSubset(ind, sub, akt.ind)
            split <- splitSubs(newSubs, indexvars, akt.ind)
            structureNeu <- sumSplitSubs(split, structure, str)
            return(structureNeu)
        }

        # create full table (every possible combination of hierarchical variables)
        structure <- expand.grid(dims)
        structure$val <- -1
        structure$index <- 0

        # we fill the compete dataset with values from "sub"
        structure <- fillStructure(sub, structure, indexvars)

        level <- list()
        # create levels
        for (i in 1:length(indexvars)) {
            level[[i]] <- factor(sort(as.character(unique(sub[,indexvars[i]]))))
        }

        # we need to fill the dataset "structure"
        st <- structure
        for (z in 1:length(indexvars)) {
            akt.var <- indexvars[z]
            # which levels do we need to calculate?
            dim <- dims[[akt.var]]
            lev <- level[[akt.var]]
            ll <- l[[akt.var]]

            todo <- dim[which(!dim %in% lev)]
            l_todo <- length(todo)

            for (i in l_todo:1) {
                subss <- st[which(st$val!=-1),]
                target <- as.character(todo[i])
                st <- fillVals(target, ll, subss, indexvars, z, st)
            }
        }
        ind <- which(st$val == -1)
        if(length(ind) > 1) {
            st <- st[-ind,]
        }
        return(st)
    }

	complete <- list()
    uniques <- list()
    dims <- list()
    toadd <- list()
    z <- 1
    for (i in indexvars) {
        minimalData[,i] <-  as.character(minimalData[,i])
        uniques[[z]] <- as.character(sort(unique(minimalData[,i])))
        complete <- genCompleteHier(minimalData, uniques[[z]], l[[z]], 1)
        sub <- complete$sub
        dims[[z]] <- complete$completeHier
        toadd <- unlist(complete$hierToAdd)
        rm(complete)
        z <- z + 1
    }

	# we need to set the colname "val" for the corresponding column in order that the functions works properly(quick'n dirty temporary solution)
	valCol <- which(! 1:ncol(sub) %in% indexvars)
	colnames(sub)[valCol] <- "val"

    # create the entire complete dataset
    fd <- createFullHier(sub, dims, l, indexvars)
    colnames(fd)[ncol(fd)] <- "geh"
    fd[,"geh"] <- ""

    if(suppVals==TRUE) {
        if(is.null(suppLimit) | is.null(suppZeros)) {
            stop("both parameters suppLimit and suppZeros are necessary!")
        }
        else {
            ### create primary suppressions
            fd[which(fd$val <= suppLimit), "geh"] <- "P"
            if(suppZeros == FALSE) {
                ind <- which(fd$val == 0)
                if(length(ind) > 0) {
                    fd[which(fd$val == 0), "geh"] <- ""
                }
            }
        }
    }

    #fd <- recodeIndexVars(fd, indexvars)
    fullData <- list()
    fullData$data <- fd
    fullData$dims <- l
    fullData$indexvars <- indexvars
    fullData$supps2check <- which(fd$geh != "")
    fullData$numberindexvars <- length(indexvars)
    class(fullData) <- "fullData"
    return(fullData)
}
