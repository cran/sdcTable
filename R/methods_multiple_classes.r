#' @aliases calc.multiple,character,list-method
#' @rdname calc.multiple-method
setMethod(f = "calc.multiple",
  signature = c("character", "list"),
  definition = function(type, input) {
    .SD <- ID <- NULL
    if (!type %in% c("makePartitions",
                     "makeAttackerProblem",
                     "calcFullProblem")) {
      stop("calc.multiple:: argument 'type' is not valid!\n")
    }

    if (type == "makePartitions") {
      return(c_make_partitions(input))
    }
    if (type == "makeAttackerProblem") {
      return(c_make_att_prob(input))
    }
    if (type == "calcFullProblem") {
      return(c_calc_full_prob(input))
    }
  }
)

setMethod("c_make_partitions", signature=c("list"), definition=function(input) {
  pI <- input$objectA
  dimInfoObj <- input$objectB
  dimInfo <- g_dim_info(dimInfoObj)
  strIDs <- g_strID(pI)

  ## create classes and groups
  tmpDat <- expand.grid(lapply(1:length(dimInfo), function(x) { 1:g_nr_levels(dimInfo[[x]]) } ))
  groups <- apply(tmpDat, 1, function(x) { paste(x, collapse="-")})
  classes <- apply(tmpDat, 1, sum)
  sortOrder <- order(classes)
  classes <- classes[sortOrder]
  classesUnique <- unique(classes)
  groups <- groups[sortOrder]
  splitGroups <- split(groups, classes)

  ## create tables for all classes and groups
  final <- list()
  final$groups <- as.list(groups)
  final$indices <- list()

  # default_codes and levels
  default_codes <- lapply(1:length(dimInfo), function(x) {
    g_default_codes(dimInfo[[x]])
  })
  dim_levels <- lapply(1:length(dimInfo), function(x) {
    g_levels(dimInfo[[x]])
  })

  # data.table to merge on
  df <- data.table(N=1:length(strIDs), strIDs=strIDs)
  setkey(df, strIDs)

  for (i in 1:length(groups)) {
    final$indices[[i]] <- list()
    levs <- as.integer(unlist(sapply(groups[[i]], strsplit, "-")))

    res <- list()
    for (z in 1:length(dimInfo)) {
      res[[z]] <- list()
      index <- which(g_levels(dimInfo[[z]]) %in% c(levs[z], levs[z]-1))
      codesDefault <- default_codes[[z]][index]
      if (levs[z] == 1) {
        res[[z]] <- codesDefault
      } else {
        levOrig <- dim_levels[[z]][index]
        diffs <- c(0,diff(levOrig))
        checkInd <- which(diffs == 1)-1
        out <- data.frame(index=index, levOrig=levOrig, codesDefault=codesDefault, ind=NA)
        out$ind[checkInd] <- 1

        checkInd <- c(checkInd, length(index))
        splitVec <- rep(0, length(index))
        for ( j in 2:length(checkInd) ) {
          if ( j < length(checkInd) ) {
            splitVec[checkInd[j-1]:(checkInd[j]-1)] <- j-1
          } else {
            splitVec[checkInd[j-1]:(checkInd[j])] <- j-1
          }
        }
        spl <- split(index, splitVec)
        counter <- 1
        for (k in 1:length(spl)) {
          rowInd <- match(spl[[k]], out$index)
          tmp <- out[rowInd,]
          if ( any(tmp[,"levOrig"]==levs[z]) ) {
            tmp <- tmp[1:(max(which(tmp$levOrig==levs[z]))),]
            res[[z]][[length(res[[z]])+1]] <- sort(unique(as.character(tmp$codesDefault)))
          }
        }
      }
    }
    final$indices[[i]] <- list()
    combs <- expand.grid(lapply(1:length(res), function(x) {
      1:length(res[[x]])
    }))

    final$indices[[i]] <- list();
    length(final$indices[[i]]) <- nrow(combs)
    for (m in 1:nrow(combs)) {
      final.strIDs <- pasteStrVec(expand(lapply(1:ncol(combs), function(x) {
        res[[x]][[combs[m,x]]]
      })), ncol(combs))
      df2 <- data.table(strIDs=final.strIDs, key="strIDs")
      final$indices[[i]][[m]] <- df[df2]$N
    }
  }
  final$nrGroups <- length(groups)
  final$nrTables <- sum(sapply(1:final$nrGroups, function(x) {
    length(final$indices[[x]])
  }))
  return(final)
})

setMethod("c_make_att_prob", signature=c("list"), definition=function(input) {
  obj <- input$objectA
  x <- slot(obj, "problemInstance")
  y <- slot(obj, "dimInfo")
  nrVars <- g_nrVars(x)

  A <- create_m_matrix(obj = obj, convert = TRUE)

  ## calculating (logical) constraints for the master problem ##
  # idea: for each constraint at least 2 suppressions must
  # exist if one xi != 0! (http://www.eia.doe.gov/ices2/missing_papers.pdf)
  newCutsMaster <- init.cutList(type='empty', input=list(nrCols=nrVars))
  #xx <- lapply(1:g_nr_rows(A), function(x) {
  # cols <- g_col_ind(g_row(A, input=list(x)))
  # v <- rep(0, nrVars)
  # v[cols] <- c(1, rep(-1, length(cols)))
  # s_add_complete_constraint(newCutsMaster) <<- list(init.cutList(type='singleCut', input=list(vals=v, dir="<=", rhs=0)))
  #})
  ################################################################

  nrConstraints <- g_nr_rows(A)
  objective <- rep(0, length=2*nrVars+nrConstraints)
  z1 <- init.simpleTriplet(type='simpleTripletDiag', input=list(nrRows=nrVars, negative=FALSE))
  z2 <- init.simpleTriplet(type='simpleTripletDiag', input=list(nrRows=nrVars, negative=TRUE))
  z <- c_bind(object=z1, input=list(z2, bindRow=FALSE))
  A <- c_bind(object=z, input=list(g_transpose(A), bindRow=FALSE))
  direction <- rep("==", g_nr_rows(A))
  rhs <- rep(0, g_nr_rows(A))

  types <- rep("C", g_nr_cols(A))
  boundsLower <- list(ind=1:g_nr_cols(A), val=c(rep(0, 2*nrVars), rep(-Inf, nrConstraints)))
  boundsUpper <- list(ind=1:g_nr_cols(A), val=c(rep(Inf, 2*nrVars), rep(Inf,  nrConstraints)))

  aProb <- new("linProb",
    objective=objective,
    constraints=A,
    direction=direction,
    rhs=rhs,
    boundsLower=boundsLower,
    boundsUpper=boundsUpper,
    types=types)
  return(list(aProb=aProb, newCutsMaster=newCutsMaster))
})

setMethod(f = "c_calc_full_prob", signature = c("list"), definition = function(input) {
  .SD <- ID <- id <- NULL
  x <- input$objectA
  y <- input$objectB
  datO <- g_raw_data(x)
  dimObj <- g_dim_info(y)

  # we have to aggregate if we are dealing with microdata
  if (g_is_microdata(x)) {
    sdvars <- setdiff(colnames(datO), key(datO))
    rawData <- datO[, lapply(.SD, sum, na.rm = TRUE), by = key(datO), .SDcols = sdvars]
  } else {
    rawData <- copy(datO)
  }
  ind.dimvars <- g_dimvar_ind(x)
  ind.freq <- g_freqvar_ind(x)

  codes <- vector("list", length = length(ind.dimvars))
  for (i in 1:length(codes)) {
    codes[[i]] <- rawData[[ind.dimvars[i]]]
    cDefault <- g_default_codes(dimObj[[i]])
    cOriginal <- g_original_codes(dimObj[[i]])
    cOriginalDups <- g_dups(dimObj[[i]])
    if (all(unique(codes[[i]]) %in% c(cOriginal, cOriginalDups))) {
      codes[[i]] <- c_match_default_codes(
        object = dimObj[[i]],
        input = rawData[[ind.dimvars[i]]]
      )
      if (sum(is.na(codes[[i]])) > 0) {
        stop(
          paste0("NA values in default codes have been generated for variable ",
          shQuote(names(dimObj)[i]),".\nPlease check the definition of this hierarchy!\n"))
      }
    } else if (all(unique(codes[[i]]) %in% cDefault)) {
      # cat("no recoding necessary!\n")
    } else {
      stop("c_calc_full_prob:: recoding not possible. Check your inputs!\n")
    }
  }

  ## calculate all possible combinations within the lowest levels of dim-vars
  ## if any combinations are missing (missing.codes), we have to set them to 0 later
  strID <- as.character(pasteStrVec(unlist(codes), length(codes)))
  exDims <- pasteStrVec(unlist(codes), length(codes))
  possDims <- sort(pasteStrVec(as.character(expand(
      lapply(dimObj, function(x) {
        g_minimal_default_codes(x)
      }), vector = TRUE
    )), length(dimObj)))
  missing.codes <- setdiff(possDims, exDims)

  ## fill the table
  nrIndexvars <- length(ind.dimvars)

  allCodes <- expand(lapply(dimObj, g_default_codes), vector = FALSE)
  fullTabObj <- data.table(ID = 1:length(allCodes[[1]]))
  for (i in 1:length(allCodes)) {
    fullTabObj[, colnames(rawData)[ind.dimvars][i] := allCodes[[i]]]
  }
  setkeyv(fullTabObj, colnames(rawData)[ind.dimvars])
  fullTabObj[, ID := NULL]

  ## revert rawData codes to default codes
  for (j in seq_along(ind.dimvars)) {
    v <- c_match_default_codes(
      object = dimObj[[j]],
      input = rawData[, get(names(dimObj)[j])]
    )
    set(rawData, NULL, names(dimObj)[j], v)
  }
  setkeyv(rawData, colnames(rawData)[ind.dimvars])

  ## replace NAs in rawData by 0 (required for aggregation)
  cols <- colnames(rawData)[(length(dimObj) + 1):ncol(rawData)]
  ind.na <- vector("list", length = length(cols))
  k <- 1
  for (j in cols) {
    ind.na[[k]] <- which(is.na(rawData[[j]]))
    set(rawData, ind.na[[k]], j, 0)
    k <- k + 1
  }

  ## merge minDat to fullDat
  fullTabObj <- merge(fullTabObj, rawData, all.x = TRUE)

  ## set missing combinations of lowest levels to 0
  ## problematic are all levels that should exist, but do not exist
  ## they are filled with 0 so that we can aggregate
  dim.vars <- colnames(fullTabObj)[ind.dimvars]
  # performance improvement
  cmd <- paste0("fullTabObj[, strID := paste0(", dim.vars[1])
  if (length(dim.vars) > 1) {
    for (i in 2:length(dim.vars)) {
      cmd <- paste0(cmd, ", ", dim.vars[i])
    }
  }
  cmd <- paste0(cmd, ")]")
  eval(parse(text = cmd))
  strID <- fullTabObj$strID
  fullTabObj[, strID := NULL]

  if (length(missing.codes) > 0) {
    index <- which(strID %in% missing.codes)
    for (i in 1:length(cols)) {
      set(fullTabObj, index, cols[i], 0)
    }
  }

  ## fill up missing dimensions
  not.finished <- TRUE

  # which indexvars have any hierarchy (not just the total?)
  # these indiecs specify the dim-variables we loop over
  useInds <- which(sapply(y@dimInfo, function(x) {
    length(x@codesOriginal) > 1
  }))

  fullTabObj[, id := .I]
  cols <- (nrIndexvars + 1):(ncol(fullTabObj) - 1)
  col.names <- names(fullTabObj)[cols]
  while (not.finished) {
    for (i in useInds) {
      if (length(dim.vars) > 1) {
        setkeyv(fullTabObj, dim.vars[-i])
      } else {
        setkeyv(fullTabObj, dim.vars[1])
      }

      cur.dim <- dimObj[[i]]@dims
      for (j in length(cur.dim):1) {
        cur.levs <-  cur.dim[[j]]
        cmd <- paste0("out <- fullTabObj[", dim.vars[i], "%in% cur.levs[-1],]")
        eval(parse(text = cmd))
        if (length(dim.vars) == 1) {
          out <- out[, lapply(.SD, sum), .SDcols = col.names]
        } else {
          out <- out[, lapply(.SD, sum), .SDcols = col.names, by = key(out)]
        }
        cmd <- paste0("row.ind <- fullTabObj[",dim.vars[i],"==cur.levs[1],id]")
        eval(parse(text = cmd))
        for (z in col.names) {
          cmd <- paste0("fullTabObj[id %in% row.ind,",z,":=out[[z]]]")
          eval(parse(text = cmd))
        }
      }
    }
    if (!is.na(fullTabObj[1, ind.freq, with = FALSE])) {
      not.finished <- FALSE
    } else {
      message("nrMissings: ", sum(is.na(fullTabObj$freq)))
    }
  }
  fullTabObj[, id := NULL]

  nrV <- nrow(fullTabObj)
  f <- fullTabObj[[ind.freq]]
  strID <- apply(fullTabObj[, dim.vars, with = FALSE], 1, paste0, collapse = "")

  # performance improvement
  cmd <- paste0("fullTabObj[,strID:=paste0(", dim.vars[1])
  if (length(dim.vars) > 1) {
    for (i in 2:length(dim.vars)) {
      cmd <- paste0(cmd, ", ", dim.vars[i])
    }
  }
  cmd <- paste0(cmd, ")]")
  eval(parse(text = cmd))
  strID <- fullTabObj$strID
  fullTabObj[, strID := NULL]

  w <- numVarsList <- NULL
  w.ind <- g_weightvar_ind(x)
  if (!is.null(w.ind)) {
    w <- fullTabObj[[w.ind]]
  }
  n.ind <- g_numvar_ind(x)
  if (!is.null(n.ind)) {
    numVarsList <- list(); length(numVarsList) <- length(n.ind)
    for (n in 1:length(n.ind)) {
      numVarsList[[n]] <- fullTabObj[[n.ind[n]]]
    }
  }

  if (length(n.ind) > 0) {
    names(numVarsList) <- colnames(g_raw_data(x))[n.ind]
  }

  ## replace 0 in rawData by NA if they have been replaced earlier
  for (i in 1:length(ind.na)) {
    if (length(ind.na[[i]]) > 0) {
      set(rawData, ind.na[[i]], cols[i], NA)
    }
  }

  s_raw_data(x) <- list(datO)

  # w are the actual weights as used in the
  # optimization procedure; comput default bounds
  # that depend on these values
  .default_bounds <- function(cost_var) {
    nr_cells <- length(cost_var)
    ll <- list()
    ll$lb <- rep(0, nr_cells)
    ll$ub <- rep(1.5 * max(cost_var), nr_cells)
    ll$LPL <-  rep(1, nr_cells)
    ll$UPL <-  rep(1, nr_cells)
    ll$SPL <-  rep(0, nr_cells)
    ll
  }

  if (!is.null(w)) {
    cost_var <- w
  } else {
    cost_var <- f
  }
  bb <- .default_bounds(cost_var = cost_var)
  problemInstance <- new(
    "problemInstance",
    strID = strID,
    Freq = f,
    w = w,
    numVars = numVarsList,
    lb = bb$lb,
    ub = bb$ub,
    LPL = bb$LPL,
    UPL = bb$UPL,
    SPL = bb$SPL,
    sdcStatus = rep("s", nrV)
  )
  problemInstance@sdcStatus[problemInstance@Freq == 0] <- "z"
  partition <- c_make_partitions(
    input = list(
      objectA = problemInstance,
      objectB = y
    )
  )
  sdcProblem <- new("sdcProblem",
    dataObj = x,
    dimInfo = y,
    problemInstance = problemInstance,
    partition = partition,
    startI = 1,
    startJ = 1,
    indicesDealtWith = NULL
  )
  return(sdcProblem)
})
