#' @importFrom stats optim median sd var quantile chisq.test pf prop.test t.test

STR.MATCH.ID <- "match.id"
STR.TARGET.VARIABLE <- "target.variable"
STR.INDEX.VARIABLE <- "index.variable"
STR.MATCH.TYPE <- "match.type"
STR.SUPPLEMENTARY.VARIABLE <- "supplementary.target.variable"
STR.SAMPLE.SIZE <- "sample.size.variable"
STR.MINIMUM <- "min"
STR.MAXIMUM <- "max"
STR.MEDIAN <- "median"
STR.QUANTILE <- "quantile"
STR.MEAN <- "mean"
STR.PROPORTION <- "proportion"
STR.STANDARD.DEVIATION <- "sd"
STR.VARIANCE <- "var"
PTN.QUANTILE <- paste0(STR.QUANTILE, "(\\.\\d+)")

#' Construct a MAIC input matrix
#' 
#' From index patient level data and a set of target baseline characteristics,
#' construct the input matrix to the maic.
#' 
#' The \code{dictionary} is a data frame containing at least 4 vectors:
#' \itemize{
#' \item "match.id" - the name of the match, used to refer to it in the
#'              matching.variables list
#' \item "target.variable" - the name of the variable in the target values
#'                     list use to inform the matching. Use dependent on
#'                     type
#' \item "index.variable" - the name of the variable in the index data frame
#'                    to match on.
#' \item "match.type" - A string indicating the match type to use. The following
#'                values are accepted:
#'   \itemize{
#'   \item minimum - records with index values lower than the target variable 
#'             will be assigned 0 weight
#'   \item maximum - records with index values greater than the target variable 
#'             will be assigned 0 weight
#'   \item median - records with index values greater than the target variable will
#'            be assigned a value of 1, those lower 0. The target for matching
#'            will be a mean of 0.5
#'   \item quantile.X - Generalisation of the median code. records with index values
#'                greater than the target variable will be assigned a value of 
#'                1, those lower 0. The target for matching will be a mean of 
#'                0.X
#'   \item mean - records will match index value directly onto target value
#'   \item proportion - as mean, with index values encoded as 1 = true, 0 = false.
#'                If target proportion is exclusive (0 or 1 exactly) then
#'                excluded members of the index population shall receive no
#'                weighting.
#'   \item sd - a matching on the square of the index value on the sum of the
#'        square of the target mean and target standard deviation. The
#'        target mean is provided by the "supplementary.target.variable"
#'   \item var - a matching on the square of the index value on the sum of the
#'         square of the target mean and the variance specified by the target
#'         variable. The target mean is provided by the 
#'         "supplementary.target.variable"
#'   }
#'  }
#'  In addition, the following vector may be necessary:
#'  \itemize{
#'  \item "supplementary.target.variable" - The name of the variable in the target
#'                                    values list that provides e.g. the mean
#'                                    for sd and var matching.
#'  }
#'  and, for estimating some p-values on difference (e.g. in proportion)
#'  \itemize{
#'  \item "sample.size.variable" - The name of the variable in the target
#'                                 values list that provides the number of
#'                                 subjects in the sample. Only for reporting.
#'  }
#' It is possible to use these match types to match on other variables by 
#' pre-processing the input correctly.
#' 
#' Finally, the \code{matching.variables} is a list or character vector containing
#' \code{match.id}s to be acted upon in this MAIC.
#' 
#' @param index A matrix or data.frame containing patient-level data
#' @param target A list containing target summary data
#' @param dictionary A data frame containing the columns "match.id",
#'                   "target.variable", "index.variable" and "match.type"
#' @param matching.variables A character vector indicating the match.id to use
#' @param x Return subject level inputs?
#' @return An object of class \code{maic.input}
#' @example R/examples/maic.example.R
#' @export
createMAICInput <- function(index,
                      target,
                      dictionary,
                      matching.variables,
                      x = FALSE){
  index <- as.data.frame(index)
  dictionary <- as.data.frame(dictionary)
  
  # Initialise variables required for return object
  target.values <- list()
  excluded <- rep(FALSE, nrow(index))
  pld.inputs <- data.frame(dummy = rep(NA, nrow(index)))
  n.adjustments <- 0
  
  if (length(matching.variables) > 0){
    ##### Sanity checking #####
    # Check vital columns in dictionary
    colnames(dictionary) <- tolower(colnames(dictionary))
    if (!all(c(STR.MATCH.ID,
               STR.TARGET.VARIABLE,
               STR.INDEX.VARIABLE,
               STR.MATCH.TYPE) %in% colnames(dictionary))){
      stop (paste("Dictionary must contain variables with the names", 
                  paste(c(STR.MATCH.ID,
                          STR.TARGET.VARIABLE,
                          STR.INDEX.VARIABLE,
                          STR.MATCH.TYPE), collapse= ", ")))
    }
    rownames(dictionary) <- dictionary[, STR.MATCH.ID]
    
    if (typeof(matching.variables) == "list"){
      matching.variables <- as.character(unlist(matching.variables))
    }
    
    ##### Start collating the data #####
    for (mv in matching.variables){
      if (!(mv %in% rownames(dictionary))){
        stop(paste(mv, "not specified in dictionary"))
      }
      target.var <- dictionary[mv, STR.TARGET.VARIABLE]
      index.var <- dictionary[mv, STR.INDEX.VARIABLE]
      match.type <- dictionary[mv, STR.MATCH.TYPE]
      
      if (!(target.var %in% names(target))){
        stop(paste(target.var, "not in target row"))
      }
      
      if (!(index.var %in% colnames(index))){
        stop(paste(index.var, "not in index data"))
      }
      
      # Could move to a switch form in future
      if (match.type == STR.MINIMUM){
        i.v <- as.numeric(index[, index.var])
        t.v <- as.numeric(target[[target.var]])
        if (!is.finite(t.v)) next()
        excluded[!is.finite(i.v) | i.v < t.v] <- TRUE
        n.adjustments <- n.adjustments + 1
      } else if (match.type == STR.MAXIMUM){
        i.v <- as.numeric(index[, index.var])
        t.v <- as.numeric(target[[target.var]])
        if (!is.finite(t.v)) next()
        excluded[!is.finite(i.v) | i.v > t.v] <- TRUE
        n.adjustments <- n.adjustments + 1
      } else if (match.type == STR.MEDIAN){
        i.v <- as.numeric(index[, index.var])
        t.v <- as.numeric(target[[target.var]])
        if (!is.finite(t.v)) next()
        excluded[!is.finite(i.v)] <- TRUE
        pld.inputs[, mv] <- ifelse(i.v < t.v, 1, 0)
        target.values[[mv]] <- 0.5
        # Check if we have any values within bounds - else fitting will fail
        if (!all(excluded) && all(pld.inputs[!excluded, mv] > 0)){
          stop("Cannot match median for ", mv, ", target is above maximum")
        } else if (!all(excluded) && all(pld.inputs[!excluded, mv] < 1)){
          stop("Cannot match median for ", mv, ", target is below minimum")
        }
        n.adjustments <- n.adjustments + 1
      } else if (match.type == STR.MEAN){
        i.v <- as.numeric(index[, index.var])
        t.v <- as.numeric(target[[target.var]])
        if (!is.finite(t.v)) next()
        excluded[!is.finite(i.v)] <- TRUE
        pld.inputs[, mv] <- i.v
        target.values[[mv]] <- t.v
        n.adjustments <- n.adjustments + 1
      } else if (match.type == STR.PROPORTION){
        i.v <- as.numeric(index[, index.var])
        t.v <- as.numeric(target[[target.var]])
        if (!is.finite(t.v)) next()
        excluded[!is.finite(i.v)] <- TRUE
        if (t.v == 0){
          excluded[i.v == 1] <- TRUE
        } else if (t.v == 1){
          excluded[i.v == 0] <- TRUE
        } else {
          pld.inputs[, mv] <- i.v
          target.values[[mv]] <- t.v
        }
        n.adjustments <- n.adjustments + 1
      } else if (match.type == STR.STANDARD.DEVIATION){
        if (!STR.SUPPLEMENTARY.VARIABLE %in% colnames(dictionary)){
          stop(paste(STR.SUPPLEMENTARY.VARIABLE, "column must be present in dictionary to use", STR.STANDARD.DEVIATION))
        }
        supp.var <- dictionary[mv, STR.SUPPLEMENTARY.VARIABLE]
        i.v <- as.numeric(index[, index.var])
        t.v <- as.numeric(target[[target.var]])
        s.v <- as.numeric(target[[supp.var]])
        if (!is.finite(t.v) || ! is.finite(s.v)) next()
        excluded[!is.finite(i.v)] <- TRUE
        pld.inputs[, mv] <- i.v * i.v
        target.values[[mv]] <- s.v * s.v + t.v * t.v
        n.adjustments <- n.adjustments + 1
      } else if (match.type == STR.VARIANCE){
        if (!STR.SUPPLEMENTARY.VARIABLE %in% colnames(dictionary)){
          stop(paste(STR.SUPPLEMENTARY.VARIABLE, "column must be present in dictionary to use", STR.VARIANCE))
        }
        supp.var <- dictionary[mv, STR.SUPPLEMENTARY.VARIABLE]
        i.v <- as.numeric(index[, index.var])
        t.v <- as.numeric(target[[target.var]])
        s.v <- as.numeric(target[[supp.var]])
        if (!is.finite(t.v) || ! is.finite(s.v)) next()
        excluded[!is.finite(i.v)] <- TRUE
        pld.inputs[, mv] <- i.v * i.v
        target.values[[mv]] <- s.v * s.v + t.v
        n.adjustments <- n.adjustments + 1
      } else if (grepl(PTN.QUANTILE, match.type)){
        mtch <- regexec(PTN.QUANTILE, match.type)
        v1 <- mtch[[1]]
        d <- as.numeric(substr(match.type,v1[2],(v1[2] + attr(v1, "match.length")[2] - 1)))
        i.v <- as.numeric(index[, index.var])
        t.v <- as.numeric(target[[target.var]])
        if (!is.finite(t.v)) next()
        excluded[!is.finite(i.v)] <- TRUE
        pld.inputs[, mv] <- ifelse(
          i.v == t.v,
          0.5,
          ifelse(i.v < t.v, 1, 0))
        target.values[[mv]] <- d
        # Check if we have any values within bounds - else fitting will fail
        if (!all(excluded) && (target.values[[mv]] > 0) && all(pld.inputs[!excluded, mv] > 0)){
          stop("Cannot match quantile for ", mv, ", target is above maximum")
        } else if (!all(excluded) && (target.values[[mv]] < 1) && all(pld.inputs[!excluded, mv] < 1)){
          stop("Cannot match quantile for ", mv, ", target is below minimum")
        }
        n.adjustments <- n.adjustments + 1
      } else {
        stop (paste(match.type, "is an unrecognised match type"))
      }
    }
    
    pld.inputs <- pld.inputs[, -which(names(pld.inputs) %in% c("dummy")), drop = FALSE]
  }
  
  
  n.matches <- length (target.values)
  n.excluded <- sum(excluded)
  
  res <- list()
  res[["n.adjustments"]] <- n.adjustments
  res[["n.matches"]] <- n.matches
  res[["excluded"]] <- excluded
  res[["target.values"]] <- target.values
  if (x){
    res[["x"]] <- pld.inputs
  }
  
  if (n.matches == 0 || n.excluded == nrow(pld.inputs)){
    # Can't go any further making the matrix, so finish constructing the return object
    res[["input.matrix"]] <- matrix(nrow = sum(!excluded),
                                    ncol = 0)
  } else {
    pld.inputs <- pld.inputs[!excluded, , drop = FALSE]
    
    nri <- nrow(pld.inputs)
    nms <- names(target.values)
    input.df <- NULL
    for (nm in nms){
      if (is.null(input.df)){
        input.df <- data.frame(pld.inputs[, nm] - target.values[[nm]])
        colnames(input.df) <- nm
      } else {
        input.df[, nm] <- pld.inputs[, nm] - target.values[[nm]]
      }
    }
    
    input.mt <- as.matrix(input.df)
    
    res[["input.matrix"]] <- input.mt
  }
  
  class(res) <- "MaicInput"
  
  return(res)
}

# Methods for the MaicWeights
#' @export
print.MaicWeight <- function(x, ...){
  print(as.numeric(x))
}

#' @export
plot.MaicWeight <- function(x, ...){
  CL <- as.list(match.call(expand.dots = TRUE))[-1]
  CL$plot <- TRUE
  #CL$x <- as.numeric(x)
  if (is.null(CL$main)){
    CL$main <- "Histogram of weights"
  }
  if (is.null(CL$xlab)){
    CL$xlab <- "Subject weight"
  }
  do.call(graphics::hist, CL)
}

# Generic function for maic weighting
#' Calculate MAIC weights
#' 
#' This function calculates the weights to apply to records for 
#' Matching-Adjusted Indirect Comparison (MAIC), from either a raw input
#' matrix or a \code{maic.input} object
#' 
#' @param x Either a \code{maic.input} object or a MAIC input matrix
#' @param opt return the optim object as attribute
#' @param keep.x return the input matrix as an attribute
#' @param ... Optional arguments to \code{\link{optim}}
#' @return A numeric vector of weights corresponding to the rows in the input
#'         matrix
#' @example R/examples/maic.example.R
#' @export
maicWeight <- function(x,
                       opt = TRUE,
                       keep.x = TRUE,
                       ...){
  UseMethod("maicWeight", x)
}

#' @export
maicWeight.MaicInput <- function(x,
                                 opt = TRUE,
                                 keep.x = TRUE,
                                 ...){
  if (x$n.matches == 0){
    z <- ifelse(x$excluded, 0, 1)
    if (opt) attr(z, "opt") <- NULL
    if (keep.x){
      attr(z, "x") <- x
    }
    attr(z, "ESS") <- sum(z)
    class(z) <- "MaicWeight"
    return(ifelse(x$excluded, 0, 1))
  }
  z <- rep(0, length(x$excluded))
  if (keep.x){
    attr(z, "x") <- x
  }
  class(z) <- "MaicWeight"
  if (all(x$excluded)){
    if (opt){
      attr(z, "opt") <- NULL
    }
    attr(z, "ESS") <- 0
    return(z)
  }
  wts <- maicWeight.default(x[["input.matrix"]],
                            opt = opt,
                            keep.x = FALSE,
                            ...)
  z[!x$excluded] <- as.numeric(wts)
  if (opt){
    attr(z, "opt") <- attr(wts, "opt")
  }
  attr(z, "ESS") <- sum(z)^2 / sum(z^2)
  z
}

#' @export
maicWeight.default <- function(x,
                               opt = TRUE,
                               keep.x = TRUE,
                               ...){
  # First, let us determine if it is going to be possible
  # to optimise on this x
  if (nrow(x) == 0){
    wt <- numeric(0)
    if (opt){
      attr(wt, "opt") <- NULL
    }
    class(wt) <- "MaicWeight"
    return(wt)
  }
  
  if (ncol(x) == 0){
    wt <- rep(1, nrow(x))
    if (opt){
      attr(wt, "opt") <- NULL
    }
    class(wt) <- "MaicWeight"
    return(wt)
  }
  
  for (cidx in seq_len(ncol(x))){
    if (all(x[, cidx] < 0) ||
        all(x[, cidx] > 0)){
      wt <- rep(0, nrow(x))
      if (opt){
        attr(wt, "opt") <- NULL
      }
      class(wt) <- "MaicWeight"
      return(wt)
    }
  }
  
  # The maic functions, as per NICE DSU
  # Objective function
  objfn <- function(a1, X){
    sum(exp(X %*% a1))
  }
  
  # Gradient function
  gradfn <- function(a1,X){
    colSums(sweep(X, 1, exp(X %*% a1), "*"))
  }
  
  optim.args <- list(...)
  if (is.null(optim.args$method)){
    optim.args$method <- "BFGS"
  }
  
  optim.args <- c(optim.args,
                  list(par = rep(0, ncol(x)),
                       fn = objfn,
                       gr = gradfn,
                       X = x))
  
  opt1 <- do.call("optim", 
                  optim.args)
  
  if (opt1$convergence != 0){
    warning("Optimisation has not converged, optim convergence code ", opt1$convergence,
            ". Please ensure you have sufficient overlapping observations",
            "for the number of matches")
  }
  
  a1 <- opt1$par
  wt <- exp(x %*% a1)
  
  if (opt){
    attr(wt, "opt") <- opt1
  }
  
  attr(wt, "ESS") <- sum(wt)^2 / sum(wt^2)
  class(wt) <- "MaicWeight"
  return (wt)
}

# Report the rebalanced covariates
#' Calculate the rebalanced covariates
#' 
#' This function calculates the raw, target and achieved covariates given
#' a set of weights.
#' Note that for mean values, bootstrapped standard errors are used and so
#' downstream values (such as p-values for difference) may differ from run
#' to run if the random number stream is not consistent
#' 
#' @param index A matrix or data.frame containing patient-level data
#' @param target A list containing target summary data
#' @param dictionary A data frame containing the columns "match.id",
#'                   "target.variable", "index.variable" and "match.type"
#' @param matching.variables A character vector indicating the match.id to use
#' @param weights A numeric vector with weights corresponding to the index 
#'                data rows
#' @param tidy A boolean - return as a data frame (otherwise list)
#' @param var.method Estimator type passed through to \code{\link{wtd.var}}.
#'                   Defaults to \code{ML}, as Bessel's correction not used in
#'                   weights generation.
#' @return An object of class \code{maic.covariates}
#' @example R/examples/maic.example.R
#' @export
reportCovariates <- function(index,
                             target,
                             dictionary,
                             matching.variables,
                             weights,
                             tidy = TRUE,
                             var.method = c("ML", "unbiased")){
  index <- as.data.frame(index)
  dictionary <- as.data.frame(dictionary)
  
  ##### Sanity checking #####
  # Check vital columns in dictionary
  colnames(dictionary) <- tolower(colnames(dictionary))
  if (!all(c(STR.MATCH.ID,
             STR.TARGET.VARIABLE,
             STR.INDEX.VARIABLE,
             STR.MATCH.TYPE) %in% colnames(dictionary))){
    stop (paste("Dictionary must contain variables with the names", 
                paste(c(STR.MATCH.ID,
                        STR.TARGET.VARIABLE,
                        STR.INDEX.VARIABLE,
                        STR.MATCH.TYPE), collapse= ", ")))
  }
  rownames(dictionary) <- dictionary[, STR.MATCH.ID]
  
  if (typeof(matching.variables) == "list"){
    matching.variables <- as.character(unlist(matching.variables))
  }
  
  ##### Start collating the data #####
  raw.value <- list()
  target.value <- list()
  adjusted.value <- list()
  unweighted.p.value <- list()
  weighted.p.value <- list()
  
  ess <- sum(weights)^2 / sum(weights^2)
  weights <- (weights/sum(weights)) * ess
  
  for (mv in matching.variables){
    if (!(mv %in% rownames(dictionary))){
      stop(paste(mv, "not specified in dictionary"))
    }
    target.var <- dictionary[mv, STR.TARGET.VARIABLE]
    index.var <- dictionary[mv, STR.INDEX.VARIABLE]
    match.type <- dictionary[mv, STR.MATCH.TYPE]
    
    if (!(target.var %in% names(target))){
      stop(paste(target.var, "not in target row"))
    }
    
    if (!(index.var %in% colnames(index))){
      stop(paste(index.var, "not in index data"))
    }
    
    # Could do a switch here, but I don't trust side effects in R
    if (match.type == STR.MINIMUM){
      i.v <- as.numeric(index[, index.var])
      t.v <- as.numeric(target[[target.var]])
      raw.value[mv] <- min(i.v, na.rm=TRUE)
      target.value[mv] <- t.v
      adjusted.value[mv] <- min(i.v * ifelse(weights > 0, 1, NA), na.rm = TRUE)
      unweighted.p.value[mv] <- NA
      weighted.p.value[mv] <- NA
    } else if (match.type == STR.MAXIMUM){
      i.v <- as.numeric(index[, index.var])
      t.v <- as.numeric(target[[target.var]])
      raw.value[mv] <- max(i.v, na.rm=TRUE)
      target.value[mv] <- t.v
      adjusted.value[mv] <- max(i.v * ifelse(weights > 0, 1, 0), na.rm = TRUE)
      unweighted.p.value[mv] <- NA
      weighted.p.value[mv] <- NA
    } else if (match.type == STR.MEDIAN){
      i.v <- as.numeric(index[, index.var])
      t.v <- as.numeric(target[[target.var]])
      raw.value[mv] <- median(i.v, na.rm = TRUE)
      target.value[mv] <- t.v
      adjusted.value[mv] <- matrixStats::weightedMedian(i.v, weights, na.rm = TRUE)
      # To test a difference in medians by the Mood's median test,
      # we need to form a combined median. We cannot do this with
      # summary data
      unweighted.p.value[mv] <- NA
      weighted.p.value[mv] <- NA
    } else if (match.type == STR.MEAN){
      i.v <- as.numeric(index[, index.var])
      t.v <- as.numeric(target[[target.var]])
      raw.value[mv] <- mean(i.v, na.rm = TRUE)
      target.value[mv] <- t.v
      adjusted.value[mv] <- Hmisc::wtd.mean(i.v, weights, na.rm = TRUE)
      
      if (is.finite(target.value[[mv]]) && is.finite(raw.value[[mv]])){
        # Caution - these are one-sample tests as cannot use
        # information about comparator sampling distribution
        unwt.tst <- t.test(i.v, mu = t.v)
        unweighted.p.value[mv] <- unwt.tst$p.value
        
        wt.tst <- weights::wtd.t.test(i.v, y = t.v, weight = weights,
                                      bootse = TRUE)
        weighted.p.value[mv] <- wt.tst$coefficients["p.value"]
      } else {
        unweighted.p.value[mv] <- NA
        weighted.p.value[mv] <- NA
      }
    } else if (match.type == STR.PROPORTION){
      i.v <- as.numeric(index[, index.var])
      t.v <- as.numeric(target[[target.var]])
      raw.value[mv] <- mean(i.v, na.rm = TRUE)
      target.value[mv] <- t.v
      adjusted.value[mv] <- sum(i.v * weights, na.rm = TRUE) / sum(weights)
      
      if (is.finite(target.value[[mv]] && is.finite(raw.value[[mv]]))){
        if (t.v <= 0 || t.v >= 1){
          unweighted.p.value[mv] <- NA
          weighted.p.value[mv] <- NA
        } else {
          if (STR.SAMPLE.SIZE %in% colnames(dictionary) &&
              !is.na(dictionary[mv, STR.SAMPLE.SIZE]) &&
              !(as.character(dictionary[mv, STR.SAMPLE.SIZE]) == "")){
            n <- as.integer(target[[as.character(dictionary[mv, STR.SAMPLE.SIZE])]])
            
            if (sum(i.v, na.rm = TRUE) < 20 || round(t.v * n) < 20 ||
                max(0, ess - sum(i.v * weights, na.rm = TRUE)) < 20 || round(n *(1-t.v)) < 20){
              sim.p <- FALSE
            } else {
              sim.p <- TRUE
            }
            
            unwt.tst <- prop.test(matrix(c(sum(i.v, na.rm = TRUE), length(i.v) - sum(i.v, na.rm = TRUE),
                                           round(t.v * n), round(n *(1-t.v))),
                                         ncol = 2, byrow = TRUE))
            unweighted.p.value[mv] <- unwt.tst$p.value
            
            wt.tst <- chisq.test(matrix(c(sum(i.v * weights, na.rm = TRUE), max(0, ess - sum(i.v * weights, na.rm = TRUE)),
                                          round(t.v * n), round(n *(1-t.v))),
                                        ncol = 2, byrow = TRUE), 
                                 simulate.p.value = sim.p)
            weighted.p.value[mv] <- wt.tst$p.value
          } else {
            # One-sample test. CAUTION!
            warning(paste("One-sample test for variable", mv))
            unwt.tst <- prop.test(sum(i.v, na.rm = TRUE), 
                                  length(i.v),
                                  t.v)
            unweighted.p.value[mv] <- unwt.tst$p.value
            
            wt.tst <- prop.test(sum(i.v * weights, na.rm = TRUE),
                                sum(weights),
                                t.v)
            weighted.p.value[mv] <- wt.tst$p.value
          }
        }
      } else {
        unweighted.p.value[mv] <- NA
        weighted.p.value[mv] <- NA
      }
    } else if (match.type == STR.STANDARD.DEVIATION){
      if (missing(var.method)) var.method <- "ML"
      i.v <- as.numeric(index[, index.var])
      t.v <- as.numeric(target[[target.var]])
      raw.value[mv] <- sd(i.v, na.rm = TRUE)
      target.value[mv] <- t.v
      adjusted.value[mv] <- sqrt(Hmisc::wtd.var(i.v, weights, na.rm = TRUE,
                                                method = var.method))
      
      if (is.finite(target.value[[mv]] && is.finite(raw.value[[mv]]))){
        # Test difference in variances using F test. Caution - requires Normality
        if (STR.SAMPLE.SIZE %in% colnames(dictionary) &&
            !is.na(dictionary[mv, STR.SAMPLE.SIZE]) &&
            !(as.character(dictionary[mv, STR.SAMPLE.SIZE]) == "")){
          n <- as.integer(target[[as.character(dictionary[mv, STR.SAMPLE.SIZE])]])
          
          ##### Unweighted test #####
          # Copy of var.test
          STATISTIC <- var(i.v, na.rm = TRUE) / (t.v ^ 2)
          PVAL <- stats::pf(STATISTIC,
                            length(i.v) - 1,
                            n - 1)
          # Use two-sided
          PVAL <- 2 * min(PVAL, 1 - PVAL)
          unweighted.p.value[mv] <- PVAL
          
          ##### Weighted test #####
          STATISTIC <- Hmisc::wtd.var(i.v, weights, na.rm = TRUE,
                                      method = var.method) / (t.v ^ 2)
          PVAL <- stats::pf(STATISTIC,
                            ess - 1,
                            n - 1)
          # Use two-sided
          PVAL <- 2 * min(PVAL, 1 - PVAL)
          weighted.p.value[mv] <- PVAL
        } else {
          # One-sample test. CAUTION!
          warning(paste("One-sample test for variable", mv))
          ##### Unweighted test #####
          # Copy of var.test
          STATISTIC <- var(i.v, na.rm = TRUE) / (t.v ^ 2)
          PVAL <- stats::pf(STATISTIC,
                            length(i.v) - 1,
                            Inf)
          # Use two-sided
          PVAL <- 2 * min(PVAL, 1 - PVAL)
          unweighted.p.value[mv] <- PVAL
          
          ##### Weighted test #####
          STATISTIC <- Hmisc::wtd.var(i.v, weights, na.rm = TRUE,
                                      method = var.method) / (t.v ^ 2)
          PVAL <- stats::pf(STATISTIC,
                            ess - 1,
                            Inf)
          # Use two-sided
          PVAL <- 2 * min(PVAL, 1 - PVAL)
          weighted.p.value[mv] <- PVAL
        }
      } else {
        unweighted.p.value[mv] <- NA
        weighted.p.value[mv] <- NA
      }
    } else if (match.type == STR.VARIANCE){
      if (missing(var.method)) var.method <- "ML"
      i.v <- as.numeric(index[, index.var])
      t.v <- as.numeric(target[[target.var]])
      raw.value[mv] <- var(i.v, na.rm = TRUE)
      target.value[mv] <- t.v
      adjusted.value[mv] <- Hmisc::wtd.var(i.v, weights, na.rm = TRUE,
                                           method = var.method)
      
      if (is.finite(target.value[[mv]] && is.finite(raw.value[[mv]]))){
        # Test difference in variances using F test. Caution - requires Normality
        if (STR.SAMPLE.SIZE %in% colnames(dictionary) &&
            !is.na(dictionary[mv, STR.SAMPLE.SIZE]) &&
            !(as.character(dictionary[mv, STR.SAMPLE.SIZE]) == "")){
          n <- as.integer(target[[as.character(dictionary[mv, STR.SAMPLE.SIZE])]])
          
          ##### Unweighted test #####
          # Copy of var.test
          STATISTIC <- var(i.v, na.rm = TRUE) / t.v
          PVAL <- stats::pf(STATISTIC,
                            length(i.v) - 1,
                            n - 1)
          # Use two-sided
          PVAL <- 2 * min(PVAL, 1 - PVAL)
          unweighted.p.value[mv] <- PVAL
          
          ##### Weighted test #####
          STATISTIC <- Hmisc::wtd.var(i.v, weights, na.rm = TRUE,
                                      method = var.method) / t.v
          PVAL <- stats::pf(STATISTIC,
                            ess - 1,
                            n - 1)
          # Use two-sided
          PVAL <- 2 * min(PVAL, 1 - PVAL)
          weighted.p.value[mv] <- PVAL
        } else {
          # One-sample test. CAUTION!
          warning(paste("One-sample test for variable", mv))
          ##### Unweighted test #####
          # Copy of var.test
          STATISTIC <- var(i.v, na.rm = TRUE) / t.v
          PVAL <- stats::pf(STATISTIC,
                            length(i.v) - 1,
                            Inf)
          # Use two-sided
          PVAL <- 2 * min(PVAL, 1 - PVAL)
          unweighted.p.value[mv] <- PVAL
          
          ##### Weighted test #####
          STATISTIC <- Hmisc::wtd.var(i.v, weights, na.rm = TRUE,
                                      method = var.method) / t.v
          PVAL <- stats::pf(STATISTIC,
                            ess - 1,
                            Inf)
          # Use two-sided
          PVAL <- 2 * min(PVAL, 1 - PVAL)
          weighted.p.value[mv] <- PVAL
        }
      } else {
        unweighted.p.value[mv] <- NA
        weighted.p.value[mv] <- NA
      }
    } else if (grepl(PTN.QUANTILE, match.type)){
      mtch <- regexec(PTN.QUANTILE, match.type)
      v1 <- mtch[[1]]
      d <- as.numeric(substr(match.type,v1[2],(v1[2] + attr(v1, "match.length")[2] - 1)))
      i.v <- as.numeric(index[, index.var])
      t.v <- as.numeric(target[[target.var]])
      raw.value[mv] <- quantile(i.v, d  / 10^(floor(log10(d)+1)), na.rm = TRUE)
      target.value[mv] <- t.v
      adjusted.value[mv] <- Hmisc::wtd.quantile(i.v, weights, probs = d  / 10^(floor(log10(d)+1)), normwt = FALSE, na.rm = TRUE)
      
      # Similar to median by Mood's test, we need to form a combined quantile. 
      # We cannot do this with summary data
      unweighted.p.value[mv] <- NA
      weighted.p.value[mv] <- NA
    } else {
      stop (paste(match.type, "is an unrecognised match type"))
    }
  }
    
  if (tidy){
    res <- data.frame(
      "matching.variable" = names(raw.value),
      "target.value" = as.numeric(target.value),
      "unadjusted.value" = as.numeric(raw.value),
      "unadjusted.delta" = as.numeric(raw.value) - as.numeric(target.value),
      "unadjusted.p.value" = as.numeric(unweighted.p.value),
      "adjusted.value" = as.numeric(adjusted.value),
      "adjusted.delta" = as.numeric(adjusted.value) - as.numeric(target.value),
      "adjusted.p.value" = as.numeric(weighted.p.value),
      stringsAsFactors = FALSE,
      row.names = names(raw.value)
    )
  } else {
    res <- list(
      "target.value" = target.value,
      "unadjusted.value" = raw.value,
      "unadjusted.delta" = as.numeric(raw.value) - as.numeric(target.value),
      "unadjusted.p.value" = unweighted.p.value,
      "adjusted.value" = adjusted.value,
      "adjusted.delta" = as.numeric(adjusted.value) - as.numeric(target.value),
      "adjusted.p.value" = weighted.p.value
    )
  }
  
  class(res) <- c(class(res), "MaicAggregates")
  
  return(res)
}

#' calculate MAIC weights
#' 
#' From index patient level data and a set of target baseline characteristics,
#' calculate MAIC weights.
#' 
#' The \code{dictionary} is a data frame containing at least 4 vectors:
#' \itemize{
#' \item "match.id" - the name of the match, used to refer to it in the
#'              matching.variables list
#' \item "target.variable" - the name of the variable in the target values
#'                     list use to inform the matching. Use dependent on
#'                     type
#' \item "index.variable" - the name of the variable in the index data frame
#'                    to match on.
#' \item "match.type" - A string indicating the match type to use. The following
#'                values are accepted:
#'   \itemize{
#'   \item minimum - records with index values lower than the target variable will
#'             be discarded
#'   \item maximum - records with index values greater than the target variable will
#'             be discarded
#'   \item median - records with index values greater than the target variable will
#'            be assigned a value of 1, those lower 0, and those equal 0.5.
#'            The target for matching will be a mean of 0.5
#'   \item quantile.X - Generalisation of the median code. records with index values
#'                greater than the target variable will be assigned a value of 
#'                1, those lower 0. The target for matching will be a mean of 
#'                0.X
#'   \item mean - records will match index value directly onto target value
#'   \item proportion - as mean, with index values encoded as 1 = true, 0 = false.
#'                If target proportion is exclusive (0 or 1 exactly) then
#'                excluded members of the index population shall receive no
#'                weighting.
#'   \item sd - a matching on the square of the index value on the sum of the
#'        square of the target mean and target standard deviation. The
#'        target mean is provided by the "supplementary.target.variable"
#'   \item var - a matching on the square of the index value on the sum of the
#'         square of the target mean and the variance specified by the target
#'         variable. The target mean is provided by the 
#'         "supplementary.target.variable"
#'   }
#'  }
#'  In addition, the following vector may be necessary:
#'  \itemize{
#'  \item "supplementary.target.variable" - The name of the variable in the target
#'                                    values list that provides e.g. the mean
#'                                    for sd and var matching.
#'  }
#' It is possible to use these match types to match on other variables, e.g.
#' variance, by pre-processing the input correctly.
#' 
#' Finally, the \code{matching.variables} is a list or character vector containing
#' \code{match.id}s to be acted upon in this MAIC.
#' 
#' @param index A matrix or data.frame containing patient-level data
#' @param target A list containing target summary data
#' @param dictionary A data frame containing the columns "match.id",
#'                   "target.variable", "index.variable" and "match.type"
#' @param matching.variables A character vector indicating the match.id to use
#' @param reporting.variables A optional character vector of matches to report
#'                            upon (defaults to \code{matching.variables})
#' @return An object of class \code{MaicAnalysis}, with components \code{weights}
#'         and \code{aggregate}, containing the weights vector and the covariate
#'         aggregate data respectively
#' @example R/examples/maic.weight.example.R
#' @export
maicMatching <- function(index,
                         target,
                         dictionary,
                         matching.variables,
                         reporting.variables = NULL){
  
  if (is.null(reporting.variables)){
    reporting.variables <- matching.variables
  }
  
  ip.mat <- createMAICInput(index,
                            target,
                            dictionary,
                            matching.variables)
  
  wts <- maicWeight(ip.mat)
    
  covars <- reportCovariates(index,
                             target,
                             dictionary,
                             reporting.variables,
                             wts)
  
  res <- list(
    "weights" = wts,
    "aggregate" = covars
  )
  
  class(res) <- "MaicAnalysis"
  
  return (res)
}