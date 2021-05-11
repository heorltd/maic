context("MAIC input")
library(maic)

test_that("A valid input matrix can be made for individual comparisons", {
  
  # Use the datasets::stackloss dataset for these tests
  data("stackloss")
  stackloss$match.conc.lt.90 <- ifelse(stackloss$Acid.Conc. < 90, 1, 0)
  
  mtch.targets <- list(
    "Air.Flow.mean" = 60,
    "Air.Flow.min" = 55,
    "Air.Flow.max" = 90,
    "Air.Flow.sd" = 9,
    "Air.Flow.var" = 81,
    "Proportion.Acid.Conc.lt.90" = 0.7,
    "Water.Temp.decile.1" = 18,
    "Water.Temp.decile.9" = 25,
    "Water.Temp.median" = 22
  )
  
  mtch.dict <- data.frame(match.id = c(
    "airflow.mean", "airflow.min", "airflow.max", "airflow.sd", 
    "airflow.var", "acidconc.prop",
    "watertemp.10tile", "watertemp.90tile", "watertemp.median"
  ),
  match.type = c(
    "mean", "min", "max", "sd",
    "var", "proportion",
    "quantile.1", "quantile.90", "median"
  ),
  target.variable = c(
    "Air.Flow.mean", "Air.Flow.min", "Air.Flow.max", "Air.Flow.sd",
    "Air.Flow.var", "Proportion.Acid.Conc.lt.90",
    "Water.Temp.decile.1", "Water.Temp.decile.9", "Water.Temp.median"
  ),
  index.variable = c(
    "Air.Flow", "Air.Flow", "Air.Flow", "Air.Flow", 
    "Air.Flow", "match.conc.lt.90",
    "Water.Temp", "Water.Temp", "Water.Temp"
  ),
  supplementary.target.variable = c(
    "", "", "", "Air.Flow.mean",
    "Air.Flow.mean", "",
    "", "", ""
  ),
  stringsAsFactors = FALSE
  )
  
  # Mean
  ipmat <- createMAICInput(index = stackloss,
                           target = mtch.targets,
                           dictionary = mtch.dict,
                           matching.variables = c("airflow.mean")
                           )
  expect_is(ipmat, "MaicInput")
  expect_equal(dim(ipmat$input.matrix), c(21, 1))
  expect_match(colnames(ipmat$input.matrix), "airflow.mean")
  expect_equal(min(ipmat$input.matrix[, 1]), min(stackloss$Air.Flow) - 
                 mtch.targets[["Air.Flow.mean"]])
  expect_equal(max(ipmat$input.matrix[, 1]), max(stackloss$Air.Flow) - 
                 mtch.targets[["Air.Flow.mean"]])
  expect_true(!(any(ipmat$excluded)))
  expect_equal(ipmat$n.adjustments, 1)
  expect_equal(ipmat$n.matches, 1)
  
  # Proportion
  ipmat <- createMAICInput(index = stackloss,
                           target = mtch.targets,
                           dictionary = mtch.dict,
                           matching.variables = c("acidconc.prop"))
  expect_is(ipmat, "MaicInput")
  expect_equal(dim(ipmat$input.matrix), c(21, 1))
  expect_match(colnames(ipmat$input.matrix), "acidconc.prop")
  expect_equal(min(ipmat$input.matrix[, 1]), 
               min(stackloss$match.conc.lt.90) - 
                 mtch.targets[["Proportion.Acid.Conc.lt.90"]])
  expect_equal(max(ipmat$input.matrix[, 1]), 
               max(stackloss$match.conc.lt.90) - 
                 mtch.targets[["Proportion.Acid.Conc.lt.90"]])
  expect_true(!(any(ipmat$excluded)))
  expect_equal(ipmat$n.adjustments, 1)
  expect_equal(ipmat$n.matches, 1)
  
  # Minimum
  ipmat <- createMAICInput(index = stackloss,
                           target = mtch.targets,
                           dictionary = mtch.dict,
                           matching.variables = c("airflow.min")
  )
  expect_is(ipmat, "MaicInput")
  expect_equal(sum(ipmat$excluded),
               sum(stackloss$Air.Flow < mtch.targets[["Air.Flow.min"]]))
  expect_equal(dim(ipmat$input.matrix), c(sum(!ipmat$excluded), 0))
  expect_equal(ipmat$n.adjustments, 1)
  expect_equal(ipmat$n.matches, 0)
  
  # Maximum
  ipmat <- createMAICInput(index = stackloss,
                           target = mtch.targets,
                           dictionary = mtch.dict,
                           matching.variables = c("airflow.max")
  )
  expect_is(ipmat, "MaicInput")
  expect_equal(sum(ipmat$excluded),
               sum(stackloss$Air.Flow > mtch.targets[["Air.Flow.max"]]))
  expect_equal(dim(ipmat$input.matrix), c(sum(!ipmat$excluded), 0))
  expect_equal(ipmat$n.adjustments, 1)
  expect_equal(ipmat$n.matches, 0)
  
  # SD
  ipmat.sd <- createMAICInput(index = stackloss,
                              target = mtch.targets,
                              dictionary = mtch.dict,
                              matching.variables = c("airflow.sd")
  )
  expect_is(ipmat.sd, "MaicInput")
  expect_equal(dim(ipmat.sd$input.matrix), c(21, 1))
  expect_match(colnames(ipmat.sd$input.matrix), "airflow.sd")
  expect_true(!(any(ipmat.sd$excluded)))
  expect_equal(ipmat.sd$n.adjustments, 1)
  expect_equal(ipmat.sd$n.matches, 1)
  expect_equal_to_reference(ipmat.sd$input.matrix,
                            "ipmat.sd.input.matrix.rds")
  
  # Var
  ipmat.var <- createMAICInput(index = stackloss,
                               target = mtch.targets,
                               dictionary = mtch.dict,
                               matching.variables = c("airflow.var")
  )
  expect_is(ipmat.var, "MaicInput")
  expect_equal(dim(ipmat.var$input.matrix), c(21, 1))
  expect_match(colnames(ipmat.var$input.matrix), "airflow.var")
  expect_true(!(any(ipmat.var$excluded)))
  expect_equal(ipmat.var$n.adjustments, 1)
  expect_equal(ipmat.var$n.matches, 1)
  expect_equal_to_reference(ipmat.var$input.matrix,
                            "ipmat.var.input.matrix.rds")
  
  expect_equal(ipmat.sd$input.matrix[,1], ipmat.var$input.matrix[,1])
  
  rm(list = c("ipmat.sd",
              "ipmat.var",
              "ipmat"))
  
  # Quantiles
  ipmat <- createMAICInput(index = stackloss,
                           target = mtch.targets,
                           dictionary = mtch.dict,
                           matching.variables = c("watertemp.10tile"))
  expect_is(ipmat, "MaicInput")
  expect_equal(dim(ipmat$input.matrix), c(21, 1))
  expect_match(colnames(ipmat$input.matrix), c("watertemp.10tile"))
  expect_true(!(any(ipmat$excluded)))
  expect_equal(ipmat$n.adjustments, 1)
  expect_equal(ipmat$n.matches, 1)
  expect_equal(sum(ipmat$input.matrix>0),
               sum(stackloss$Water.Temp < 
                     mtch.targets[["Water.Temp.decile.1"]]))
  
  ipmat <- createMAICInput(index = stackloss,
                           target = mtch.targets,
                           dictionary = mtch.dict,
                           matching.variables = c("watertemp.90tile"))
  expect_is(ipmat, "MaicInput")
  expect_equal(dim(ipmat$input.matrix), c(21, 1))
  expect_match(colnames(ipmat$input.matrix), c("watertemp.90tile"))
  expect_true(!(any(ipmat$excluded)))
  expect_equal(ipmat$n.adjustments, 1)
  expect_equal(ipmat$n.matches, 1)
  expect_equal(sum(ipmat$input.matrix>0),
               sum(stackloss$Water.Temp < 
                     mtch.targets[["Water.Temp.decile.9"]]))
  
  # Median
  ipmat <- createMAICInput(index = stackloss,
                           target = mtch.targets,
                           dictionary = mtch.dict,
                           matching.variables = c("watertemp.median"))
  expect_is(ipmat, "MaicInput")
  expect_equal(dim(ipmat$input.matrix), c(21, 1))
  expect_match(colnames(ipmat$input.matrix), c("watertemp.median"))
  expect_true(!(any(ipmat$excluded)))
  expect_equal(ipmat$n.adjustments, 1)
  expect_equal(ipmat$n.matches, 1)
  expect_equal(sum(ipmat$input.matrix>0),
               sum(stackloss$Water.Temp < 
                     mtch.targets[["Water.Temp.median"]]))
})

test_that("Sensible errors are thrown if the input constructor is fed bad input", {
  # Use the datasets::stackloss dataset for these tests
  data("stackloss")
  stackloss$match.conc.lt.90 <- ifelse(stackloss$Acid.Conc. < 90, 1, 0)
  
  mtch.targets <- list(
    "Air.Flow.mean" = 60,
    "Air.Flow.min" = 55,
    "Air.Flow.max" = 90,
    "Air.Flow.sd" = 9,
    "Air.Flow.var" = 81,
    "Proportion.Acid.Conc.lt.90" = 0.7,
    "Water.Temp.decile.1" = 18,
    "Water.Temp.decile.9" = 25,
    "Water.Temp.median" = 22
  )
  
  mtch.dict <- data.frame(match.id = c(
    "airflow.mean", "airflow.min", "airflow.max", "airflow.sd", 
    "airflow.var", "acidconc.prop",
    "watertemp.10tile", "watertemp.90tile", "watertemp.median"
  ),
  match.type = c(
    "mean", "min", "max", "sd",
    "var", "proportion",
    "quantile.1", "quantile.90", "median"
  ),
  target.variable = c(
    "Air.Flow.mean", "Air.Flow.min", "Air.Flow.max", "Air.Flow.sd",
    "Air.Flow.var", "Proportion.Acid.Conc.lt.90",
    "Water.Temp.decile.1", "Water.Temp.decile.9", "Water.Temp.median"
  ),
  index.variable = c(
    "Air.Flow", "Air.Flow", "Air.Flow", "Air.Flow", 
    "Air.Flow", "match.conc.lt.90",
    "Water.Temp", "Water.Temp", "Water.Temp"
  ),
  supplementary.target.variable = c(
    "", "", "", "Air.Flow.mean",
    "Air.Flow.mean", "",
    "", "", ""
  ),
  stringsAsFactors = FALSE
  )
  
  # Check that things work for ordinary creation
  ipmat <- createMAICInput(index = stackloss,
                           target = mtch.targets,
                           dictionary = mtch.dict,
                           matching.variables = c("airflow.mean",
                                                  "airflow.sd",
                                                  "airflow.min",
                                                  "airflow.max",
                                                  "acidconc.prop",
                                                  "watertemp.10tile",
                                                  "watertemp.90tile",
                                                  "watertemp.median")
  )
  expect_is(ipmat, "MaicInput")
  
  # Check that things work for truncated variables
  ipmat <- createMAICInput(ind = stackloss,
                           tar = mtch.targets,
                           dic = mtch.dict,
                           mat = c("airflow.mean",
                                                  "airflow.sd",
                                                  "airflow.min",
                                                  "airflow.max",
                                                  "acidconc.prop",
                                                  "watertemp.10tile",
                                                  "watertemp.90tile",
                                                  "watertemp.median")
  )
  expect_is(ipmat, "MaicInput")
  
  # Check that things work if there are no matches
  ipmat <- createMAICInput(ind = stackloss,
                           tar = mtch.targets,
                           dic = mtch.dict,
                           mat = c()
  )
  expect_is(ipmat, "MaicInput")
  expect_true(!(any(ipmat$excluded)))
  expect_equal(ipmat$n.adjustments, 0)
  expect_equal(ipmat$n.matches, 0)
  
  # Check failure if there are vital colnames missing
  d2 <- mtch.dict
  d2 <- d2[, -which(colnames(d2) %in% c("match.id")), drop = FALSE]
  expect_error(createMAICInput(ind = stackloss,
                           tar = mtch.targets,
                           dic = d2,
                           mat = c("airflow.mean",
                                   "airflow.sd",
                                   "airflow.min",
                                   "airflow.max",
                                   "acidconc.prop",
                                   "watertemp.10tile",
                                   "watertemp.90tile",
                                   "watertemp.median")))
  
  # Check easy passthrough if there are no matches
  ipmat <- createMAICInput(ind = stackloss,
                               tar = mtch.targets,
                               dic = d2,
                               mat = c())
  expect_is(ipmat, "MaicInput")
  expect_true(!(any(ipmat$excluded)))
  expect_equal(ipmat$n.adjustments, 0)
  expect_equal(ipmat$n.matches, 0)
  
  # Check failure if there are aux colnames missing
  d2 <- mtch.dict
  d2 <- d2[, -which(colnames(d2) %in% c("supplementary.target.variable")),
           drop = FALSE]
  expect_error(createMAICInput(ind = stackloss,
                               tar = mtch.targets,
                               dic = d2,
                               mat = c("airflow.mean",
                                       "airflow.sd",
                                       "airflow.min",
                                       "airflow.max",
                                       "acidconc.prop",
                                       "watertemp.10tile",
                                       "watertemp.90tile",
                                       "watertemp.median")))
  
  # Check OK if we take out the SD match
  ipmat <- createMAICInput(ind = stackloss,
                           tar = mtch.targets,
                           dic = d2,
                           mat = c("airflow.mean",
                                   "airflow.min",
                                   "airflow.max",
                                   "acidconc.prop",
                                   "watertemp.10tile",
                                   "watertemp.90tile",
                                   "watertemp.median"))
  expect_is(ipmat, "MaicInput")
  
  rm(d2)
  
  # What happens if we try to set up a quantile match outside the domain
  # of our data?
  mtch.targets.2 <- mtch.targets
  mtch.targets.2[["Water.Temp.decile.1"]] <- min(stackloss$Water.Temp) - 1
  expect_error(createMAICInput(ind = stackloss,
                               tar = mtch.targets.2,
                               dic = mtch.dict,
                               mat = c("airflow.mean",
                                   "airflow.min",
                                   "airflow.max",
                                   "acidconc.prop",
                                   "watertemp.10tile",
                                   "watertemp.90tile",
                                   "watertemp.median")),
    "Cannot match quantile for watertemp.10tile, target is below minimum")
  
  mtch.targets.2 <- mtch.targets
  mtch.targets.2[["Water.Temp.decile.9"]] <- max(stackloss$Water.Temp) + 1
  expect_error(createMAICInput(ind = stackloss,
                               tar = mtch.targets.2,
                               dic = mtch.dict,
                               mat = c("airflow.mean",
                                       "airflow.min",
                                       "airflow.max",
                                       "acidconc.prop",
                                       "watertemp.10tile",
                                       "watertemp.90tile",
                                       "watertemp.median")),
    "Cannot match quantile for watertemp.90tile, target is above maximum")
  
  # Median match outside domain
  mtch.targets.2 <- mtch.targets
  mtch.targets.2[["Water.Temp.median"]] <- min(stackloss$Water.Temp) - 1
  expect_error(createMAICInput(ind = stackloss,
                               tar = mtch.targets.2,
                               dic = mtch.dict,
                               mat = c("airflow.mean",
                                       "airflow.min",
                                       "airflow.max",
                                       "acidconc.prop",
                                       "watertemp.10tile",
                                       "watertemp.90tile",
                                       "watertemp.median")),
    "Cannot match median for watertemp.median, target is below minimum")
  
  mtch.targets.2 <- mtch.targets
  mtch.targets.2[["Water.Temp.median"]] <- max(stackloss$Water.Temp) + 1
  expect_error(createMAICInput(ind = stackloss,
                               tar = mtch.targets.2,
                               dic = mtch.dict,
                               mat = c("airflow.mean",
                                       "airflow.min",
                                       "airflow.max",
                                       "acidconc.prop",
                                       "watertemp.10tile",
                                       "watertemp.90tile",
                                       "watertemp.median")),
    "Cannot match median for watertemp.median, target is above maximum")
  
  # Data missing from index
  stackloss.2 <- stackloss[, - which(colnames(stackloss) %in% "Air.Flow")]
  expect_error(createMAICInput(index = stackloss.2,
                           target = mtch.targets,
                           dictionary = mtch.dict,
                           matching.variables = c("airflow.mean",
                                                  "airflow.sd",
                                                  "airflow.min",
                                                  "airflow.max",
                                                  "acidconc.prop",
                                                  "watertemp.10tile",
                                                  "watertemp.90tile",
                                                  "watertemp.median")),
    "Air.Flow not in index data")
  
  # Data missing from target
  mtch.targets.2 <- mtch.targets[-which(names(mtch.targets) == "Air.Flow.mean")]
  expect_error(createMAICInput(index = stackloss,
                  target = mtch.targets.2,
                  dictionary = mtch.dict,
                  matching.variables = c("airflow.mean",
                                         "airflow.sd",
                                         "airflow.min",
                                         "airflow.max",
                                         "acidconc.prop",
                                         "watertemp.10tile",
                                         "watertemp.90tile",
                                         "watertemp.median")),
    "Air.Flow.mean not in target row")
})
  