# maic.R

###############################################################################
# SCRIPT:
# Name:       maic
# Date:       02 Aug 2017
# Version:    0.0.2
# Authors:    Rob Young (robert.young@heor.co.uk)
#
# Description:
#             Functions for handling the matching-adjusted indirect comparison
#             processes. In this, patient-level data from an INDEX study is
#             weighted to match a TARGET study, represented by summary data
#
# Version History:
# - Version:  0.0.3
# - Date:     28 Sep 2017
# - Author:   Rob Young
# - Changes:  Add in standard deviation / variance matching
# -----------------------------------------------------------------------------
# - Version:  0.0.2
# - Date:     02 Aug
# - Author:   Rob Young
# - Changes:  Upper bound for proportions was set at 0, not 1!
# -----------------------------------------------------------------------------
# - Version:  0.0.1
# - Date:     14 June
# - Author:   Rob Young
# - Changes:  Original version
# -----------------------------------------------------------------------------
###############################################################################

#' @importFrom stats optim median sd var quantile

STR.MATCH.ID <- "match.id"
STR.TARGET.VARIABLE <- "target.variable"
STR.INDEX.VARIABLE <- "index.variable"
STR.MATCH.TYPE <- "match.type"
STR.SUPPLEMENTARY.VARIABLE <- "supplementary.target.variable"
STR.MINIMUM <- "min"
STR.MAXIMUM <- "max"
STR.MEDIAN <- "median"
STR.QUANTILE <- "quantile"
STR.MEAN <- "mean"
STR.PROPORTION <- "proportion"
STR.STANDARD.DEVIATION <- "sd"
STR.VARIANCE <- "var"
PTN.QUANTILE <- paste0(STR.QUANTILE, "\\.(d+)")

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
#'   \item minimum - records with index values lower than the target variable will
#'             be discarded
#'   \item maximum - records with index values greater than the target variable will
#'             be discarded
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
#' @return An object of class \code{maic.input}
#' @example R/maic.example.R
#' @export
createMAICInput <- function(index,
                      target,
                      dictionary,
                      matching.variables){
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
  target.values <- list()
  excluded <- rep(FALSE, nrow(index))
  pld.inputs <- data.frame(dummy = rep(NA, nrow(index)))
  n.adjustments <- 0
  
  
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
      if (!is.finite(t.v)) next()
      excluded[i.v < t.v] <- TRUE
      n.adjustments <- n.adjustments + 1
    } else if (match.type == STR.MAXIMUM){
      i.v <- as.numeric(index[, index.var])
      t.v <- as.numeric(target[[target.var]])
      if (!is.finite(t.v)) next()
      excluded[i.v > t.v] <- TRUE
      n.adjustments <- n.adjustments + 1
    } else if (match.type == STR.MEDIAN){
      i.v <- as.numeric(index[, index.var])
      t.v <- as.numeric(target[[target.var]])
      if (!is.finite(t.v)) next()
      pld.inputs[, mv] <- ifelse(i.v < t.v, 1, 0)
      target.values[[mv]] <- 0.5
      n.adjustments <- n.adjustments + 1
    } else if (match.type == STR.MEAN){
      i.v <- as.numeric(index[, index.var])
      t.v <- as.numeric(target[[target.var]])
      if (!is.finite(t.v)) next()
      pld.inputs[, mv] <- i.v
      target.values[[mv]] <- t.v
      n.adjustments <- n.adjustments + 1
    } else if (match.type == STR.PROPORTION){
      i.v <- as.numeric(index[, index.var])
      t.v <- as.numeric(target[[target.var]])
      if (!is.finite(t.v)) next()
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
      pld.inputs[, mv] <- i.v * i.v
      target.values[[mv]] <- s.v * s.v + t.v
      n.adjustments <- n.adjustments + 1
    } else if (grepl(PTN.QUANTILE, match.type)){
      mtch <- regexec(PTN.QUANTILE, match.type)
      v1 <- mtch[[1]][1]
      d <- as.numeric(substr(match.type,v1[2],(v1[2] + attr(v1, "match.length")[2] - 1)))
      i.v <- as.numeric(index[, index.var])
      t.v <- as.numeric(target[[target.var]])
      if (!is.finite(t.v)) next()
      pld.inputs[, mv] <- ifelse(i.v < t.v, 1, 0)
      target.values[[mv]] <- d
      n.adjustments <- n.adjustments + 1
    } else {
      stop (paste(match.type, "is an unrecognised match type"))
    }
  }
  
  pld.inputs <- pld.inputs[, -which(names(pld.inputs) %in% c("dummy")), drop = FALSE]
  
  n.matches <- length (target.values)
  n.excluded <- sum(excluded)
  
  res <- list()
  res[["n.adjustments"]] <- n.adjustments
  res[["n.matches"]] <- n.matches
  res[["excluded"]] <- excluded
  
  if (n.matches == 0 || n.excluded == nrow(pld.inputs)){
    # Can't go any further making the matrix, so finish constructing the return object
    res[["input.matrix"]] <- matrix()
  } else {
    pld.inputs <- pld.inputs[!excluded, , drop = FALSE]
    
    nri <- nrow(pld.inputs)
    nms <- names(target.values)
    input.df <- NULL
    for (nm in nms){
      if (is.null(input.df)){
        input.df <- data.frame(pld.inputs[, nm] - target.values[[nm]])
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

# Constructor function for maic.input
#' Constructor for a maic.input object
#' 
#' @param n.adjustments Numeric, number of variables that have had any
#'                      comparison performed
#' @param n.matches Numeric, number of matching variables
#' @param excluded Logical vector; which index rows have been excluded from 
#'                 matching
#' @param input.matrix Numeric matrix, centred MAIC input matrix
#' @return An object of class \code{maic.input}
#' @export
MaicInput <- function(n.adjustments,
                       n.matches,
                       excluded,
                       input.matrix){
  structure(list("n.adjustments" = n.adjustments,
                 "n.matches" = n.matches,
                 "excluded" = excluded,
                 "input.matrix" = input.matrix),
            class = "maic.input")
}
  

# Generic function for maic weighting
#' Calculate MAIC weights
#' 
#' This function calculates the weights to apply to records for 
#' Matching-Adjusted Indirect Comparison (MAIC), from either a raw input
#' matrix or a \code{maic.input} object
#' 
#' @param x Either a \code{maic.input} object or a MAIC input matrix
#' @return A numeric vector of weights corresponding to the rows in the input
#'         matrix
#' @example R/maic.example.R
#' @export
maicWeight <- function(x){
  UseMethod("maicWeight", x)
}

#' Calculate MAIC weights
#' 
#' This function calculates the weights to apply to records for 
#' Matching-Adjusted Indirect Comparison (MAIC), from a \code{maic.input} 
#' object
#' 
#' @param x A \code{maic.input} object
#' @return A numeric vector of weights corresponding to the rows in the input
#'         matrix
#' @example R/maic.example.R
#' @export
maicWeight.MaicInput <- function(x){
  maicWeight.default(x[["input.matrix"]])
}

# The basic MAIC weighting code
#' Calculate MAIC weights
#' 
#' This function calculates the weights to apply to records for 
#' Matching-Adjusted Indirect Comparison (MAIC), from a raw input matrix
#' 
#' @param x A MAIC input matrix
#' @return A numeric vector of weights corresponding to the rows in the input
#'         matrix
#' @example R/maic.example.R
#' @export
maicWeight.default <- function(x){
  # The maic functions, as per NICE DSU
  # Objective function
  objfn <- function(a1, X){
    sum(exp(X %*% a1))
  }
  
  # Gradient function
  gradfn <- function(a1,X){
    colSums(sweep(X,1,exp(X %*% a1), "*"))
  }
  
  opt1 <- optim(par=rep(0, ncol(x)),
                fn=objfn,
                gr=gradfn,
                X=x,
                method="BFGS")
  a1 <- opt1$par
  wt <- exp(x %*% a1)
  
  return (wt)
}

# Report the rebalanced covariates
#' Calculate the rebalanced covariates
#' 
#' This function calculates the raw, target and achieved covariates given
#' a set of weights
#' 
#' @param index A matrix or data.frame containing patient-level data
#' @param target A list containing target summary data
#' @param dictionary A data frame containing the columns "match.id",
#'                   "target.variable", "index.variable" and "match.type"
#' @param matching.variables A character vector indicating the match.id to use
#' @param weights A numeric vector with weights corresponding to the index 
#'                data rows
#' @return An object of class \code{maic.covariates}
#' @example R/maic.example.R
#' @export
reportCovariates <- function(index,
                             target,
                             dictionary,
                             matching.variables,
                             weights){
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
      raw.value[[mv]] <- min(i.v, na.rm=TRUE)
      target.value[[mv]] <- t.v
      adjusted.value[[mv]] <- min(i.v * ifelse(weights > 0, 1, NA), na.rm = TRUE)
    } else if (match.type == STR.MAXIMUM){
      i.v <- as.numeric(index[, index.var])
      t.v <- as.numeric(target[[target.var]])
      raw.value[[mv]] <- max(i.v, na.rm=TRUE)
      target.value[[mv]] <- t.v
      adjusted.value[[mv]] <- max(i.v * ifelse(weights > 0, 1, 0))
    } else if (match.type == STR.MEDIAN){
      i.v <- as.numeric(index[, index.var])
      t.v <- as.numeric(target[[target.var]])
      raw.value[[mv]] <- median(i.v, na.rm = TRUE)
      target.value[[mv]] <- t.v
      adjusted.value[[mv]] <- matrixStats::weightedMedian(i.v, weights)
    } else if (match.type == STR.MEAN){
      i.v <- as.numeric(index[, index.var])
      t.v <- as.numeric(target[[target.var]])
      raw.value[[mv]] <- mean(i.v, na.rm = TRUE)
      target.value[[mv]] <- t.v
      adjusted.value[[mv]] <- sum(i.v * weights) / sum(weights)
    } else if (match.type == STR.PROPORTION){
      i.v <- as.numeric(index[, index.var])
      t.v <- as.numeric(target[[target.var]])
      raw.value[[mv]] <- mean(i.v, na.rm = TRUE)
      target.value[[mv]] <- t.v
      adjusted.value[[mv]] <- sum(i.v * weights) / sum(weights)
    } else if (match.type == STR.STANDARD.DEVIATION){
      i.v <- as.numeric(index[, index.var])
      t.v <- as.numeric(target[[target.var]])
      raw.value[[mv]] <- sd(i.v, na.rm = TRUE)
      target.value[[mv]] <- t.v
      adjusted.value[[mv]] <- sqrt(Hmisc::wtd.var(i.v, na.rm = TRUE))
    } else if (match.type == STR.VARIANCE){
      i.v <- as.numeric(index[, index.var])
      t.v <- as.numeric(target[[target.var]])
      raw.value[[mv]] <- var(i.v, na.rm = TRUE)
      target.value[[mv]] <- t.v
      adjusted.value[[mv]] <- Hmisc::wtd.var(i.v, na.rm = TRUE)
    } else if (grepl(PTN.QUANTILE, match.type)){
      mtch <- regexec(PTN.QUANTILE, match.type)
      v1 <- mtch[[1]][1]
      d <- as.numeric(substr(match.type,v1[2],(v1[2] + attr(v1, "match.length")[2] - 1)))
      i.v <- as.numeric(index[, index.var])
      t.v <- as.numeric(target[[target.var]])
      raw.value[[mv]] <- quantile(i.v, d, na.rm = TRUE)
      target.value[[mv]] <- t.v
      adjusted.value[[mv]] <- Hmisc::wtd.quantile(i.v, weights, normwt = TRUE, na.rm = TRUE)
    } else {
      stop (paste(match.type, "is an unrecognised match type"))
    }
  }
    
  res <- list(
    "raw.values" = raw.value,
    "target.values" = target.value,
    "adjusted.values" = adjusted.value
  )
  
  class(res) <- "MaicCovariates"
  
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
#' @param check.residuals Logical - calculate residuals to check
#' @param residual.warning.level Numeric - level at which to raise a warning
#'                               that matching has not succeeded
#' @return An object of class \code{MAICweights}
#' @example R/maic.weight.example.R
#' @export
maicMatching <- function(index,
                         target,
                         dictionary,
                         matching.variables,
                         reporting.variables = NULL,
                         check.residuals = TRUE,
                         residual.warning.level = 1e-3){
  
  if (is.null(reporting.variables)){
    reporting.variables <- matching.variables
  }
  
  ip.mat <- createMAICInput(index,
                            target,
                            dictionary,
                            matching.variables)
  
  
  if (ip.mat$n.matches == 0){
    wts <- rep(1, nrow(index))
    wts[ipmat$excluded] <- 0
  } else {
    wt <- maicWeight(ip.mat)
    wts <- rep(0, nrow(index))
    wts[!ipmat$excluded] <- wt
  }
  
  covars <- reportCovariates(index,
                             target,
                             dictionary,
                             reporting.variables,
                             wts)
  
  res <- list(
    "weights" = wts,
    "covariates" = covars
  )
  
  class(res) <- "MAICWeights"
  
  if (check.residuals){
    mtch.covars <- reportCovariates(index,
                                    target,
                                    dictionary,
                                    matching.variables,
                                    wts)
    
    resids <- numeric(0)
    
    for (mv in matching.variables){
      rownames(dictionary) <- dictionary[, STR.MATCH.ID]
      mv.type <- dictionary[mv, STR.MATCH.TYPE]
      
      resids[mv] <- 0
      
      if (!(mv.type %in% c("min", "max"))){
        tgt <- mtch.covars$target.values[[mv]]
        if (is.finite(tgt)){
          adj <- mtch.covars$adjusted.values[[mv]]
          
          resids[mv] <- (tgt - adj) ^ 2
        }
      }
    }
    
    res[["residuals"]] <- resids
    
    if (sum(resids) > residual.warning.level){
      warning("Matching failure")
    }
  }
  
  return (res)
}