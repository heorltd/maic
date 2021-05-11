context("MAIC weighting")
library(maic)

test_that("Normal processing of correctly specified matching",{
  
  data("stackloss")
  stackloss$match.conc.lt.90 <- ifelse(stackloss$Acid.Conc. < 90, 1, 0)
  
  mtch.targets <- list(
    "Air.Flow.mean" = 65,
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
  
  ##### Each matching type #####
  # Mean
  ipmat <- createMAICInput(index = stackloss,
                           target = mtch.targets,
                           dictionary = mtch.dict,
                           matching.variables = c("airflow.mean")
  )
  
  wts <- maicWeight(ipmat, control = list(reltol = 1e-12))
  expect_equal(length(wts), 21)
  
  # We expect all unique values of Air.Flow to have the same weight
  af.fct <- factor(stackloss$Air.Flow)
  same.weight <- TRUE
  for(lvl in levels(af.fct)){
    same.weight <- same.weight || length(unique(wts[af.fct==lvl])) != 1
  }
  expect_true(same.weight)
  # This is a little bit strong maybe, but this is the target:
  expect_equal(sum(wts * stackloss$Air.Flow) / sum(wts), 
               mtch.targets[["Air.Flow.mean"]],
               tolerance = 1e-6)
  
  # Proportion
  ipmat <- createMAICInput(index = stackloss,
                           target = mtch.targets,
                           dictionary = mtch.dict,
                           matching.variables = c("acidconc.prop"))
  wts <- maicWeight(ipmat)
  expect_equal(length(wts),21)
  
  # We expect all unique values of Air.Flow to have the same weight
  af.fct <- factor(stackloss$match.conc.lt.90)
  same.weight <- TRUE
  for(lvl in levels(af.fct)){
    same.weight <- same.weight || length(unique(wts[af.fct==lvl])) != 1
  }
  expect_true(same.weight)
  expect_equal(sum(wts * stackloss$match.conc.lt.90) / sum(wts), 
               mtch.targets[["Proportion.Acid.Conc.lt.90"]],
               tolerance = 1e-6)
  
  # Minimum
  ipmat <- createMAICInput(index = stackloss,
                           target = mtch.targets,
                           dictionary = mtch.dict,
                           matching.variables = c("airflow.min")
  )
  wts <- maicWeight(ipmat)
  expect_equal(length(wts), 21)
  expect_true(all(wts[stackloss$Air.Flow < mtch.targets[["Air.Flow.min"]]] == 0))
  expect_true(all(wts[stackloss$Air.Flow >= mtch.targets[["Air.Flow.min"]]] == 1))
  
  # Maximum
  ipmat <- createMAICInput(index = stackloss,
                           target = mtch.targets,
                           dictionary = mtch.dict,
                           matching.variables = c("airflow.max")
  )
  wts <- maicWeight(ipmat)
  expect_equal(length(wts), 21)
  expect_true(all(wts[stackloss$Air.Flow > 
                        mtch.targets[["Air.Flow.max"]]] == 0))
  expect_true(all(wts[stackloss$Air.Flow <= 
                        mtch.targets[["Air.Flow.max"]]] == 1))
  
  # SD
  ipmat.sd <- createMAICInput(index = stackloss,
                              target = mtch.targets,
                              dictionary = mtch.dict,
                              matching.variables = c("airflow.mean", 
                                                     "airflow.sd")
  )
  wts.sd <- maicWeight(ipmat.sd,
                    control = list(reltol = 1e-12))
  expect_equal(length(wts.sd), 21)
  
  # We expect all unique values of Air.Flow to have the same weight
  af.fct <- factor(stackloss$Air.Flow)
  same.weight <- TRUE
  for(lvl in levels(af.fct)){
    same.weight <- same.weight || length(unique(wts.sd[af.fct==lvl])) != 1
  }
  expect_true(same.weight)
  expect_equal(sqrt(Hmisc::wtd.var(stackloss$Air.Flow,
                                   wts.sd,
                                   method = "ML")),
               mtch.targets[["Air.Flow.sd"]],
               tolerance = 1e-6)
  
  # Variance
  ipmat.var <- createMAICInput(index = stackloss,
                               target = mtch.targets,
                               dictionary = mtch.dict,
                               matching.variables = c("airflow.mean",
                                                      "airflow.var"))
  wts.var <- maicWeight(ipmat.var,
                    control = list(reltol = 1e-12))
  expect_equal(length(wts.var), 21)
  
  # We expect all unique values of Air.Flow to have the same weight
  af.fct <- factor(stackloss$Air.Flow)
  same.weight <- TRUE
  for(lvl in levels(af.fct)){
    same.weight <- same.weight || length(unique(wts.var[af.fct==lvl])) != 1
  }
  expect_true(same.weight)
  expect_equal(Hmisc::wtd.var(stackloss$Air.Flow,
                                   wts.var,
                                   method = "ML"),
               mtch.targets[["Air.Flow.var"]],
               tolerance = 1e-6)
  
  expect_true(all(wts.sd == wts.var))
  
  rm(ipmat, ipmat.sd, ipmat.var, wts, wts.sd, wts.var)
  
  # Quantiles
  ipmat <- createMAICInput(index = stackloss,
                           target = mtch.targets,
                           dictionary = mtch.dict,
                           matching.variables = "watertemp.10tile")
  wts <- maicWeight(ipmat,
                    control = list(reltol = 1e-12))
  
  expect_equal(length(wts), 21)
  expect_true(length(unique(wts[stackloss$Water.Temp < 
                           mtch.targets[["Water.Temp.decile.1"]]])) == 1)
  expect_true(length(unique(wts[stackloss$Water.Temp >= 
                           mtch.targets[["Water.Temp.decile.1"]]])) == 1)
  # Have to have quite a high tolerance here
  expect_equal(as.numeric(Hmisc::wtd.quantile(stackloss$Water.Temp,
                                   wts,
                                   probs = 0.1)),
               mtch.targets[["Water.Temp.decile.1"]],
               tolerance = 1e-2)
  
  ipmat <- createMAICInput(index = stackloss,
                           target = mtch.targets,
                           dictionary = mtch.dict,
                           matching.variables = "watertemp.90tile")
  wts <- maicWeight(ipmat,
                    control = list(reltol = 1e-12))
  
  expect_equal(length(wts), 21)
  expect_true(length(unique(wts[stackloss$Water.Temp <
                           mtch.targets[["Water.Temp.decile.9"]]])) == 1)
  expect_true(length(unique(wts[stackloss$Water.Temp >= 
                           mtch.targets[["Water.Temp.decile.9"]]])) == 1)
  expect_equal(as.numeric(Hmisc::wtd.quantile(stackloss$Water.Temp,
                                              wts,
                                              probs = 0.9)),
               mtch.targets[["Water.Temp.decile.9"]],
               tolerance = 1e-2)
})

test_that("Highly stressed matching",{
  
  data("stackloss")
  stackloss$match.conc.lt.90 <- ifelse(stackloss$Acid.Conc. < 90, 1, 0)
  
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
  
  ##### Subgrouping #####
  mtch.targets <- list(
    "Proportion.Acid.Conc.lt.90" = 0
  )
  ipmat <- createMAICInput(index = stackloss,
                           target = mtch.targets,
                           dictionary = mtch.dict,
                           matching.variables = c("acidconc.prop")
  )
  wts <- maicWeight(ipmat, control = list(reltol = 1e-12))
  expect_equal(sum(wts!=0), sum(stackloss$match.conc.lt.90 == 0))
  expect_true(all(wts[stackloss$match.conc.lt.90 == 1] == 0))
  
  mtch.targets <- list(
    "Proportion.Acid.Conc.lt.90" = 1
  )
  ipmat <- createMAICInput(index = stackloss,
                           target = mtch.targets,
                           dictionary = mtch.dict,
                           matching.variables = c("acidconc.prop")
  )
  wts <- maicWeight(ipmat, control = list(reltol = 1e-12))
  expect_equal(sum(wts!=0), sum(stackloss$match.conc.lt.90 == 1))
  expect_true(all(wts[stackloss$match.conc.lt.90 == 0] == 0))
  
  ##### Force errors #####
  # Mean exceeding domain
  mtch.targets <- list(
    "Air.Flow.mean" = 85
  )
  ipmat <- createMAICInput(index = stackloss,
                           target = mtch.targets,
                           dictionary = mtch.dict,
                           matching.variables = c("airflow.mean")
  )
  
  wts <- maicWeight(ipmat, control = list(reltol = 1e-12))
  expect_equal(length(wts), 21)
  expect_true(all(wts == 0))
  
  # More targets than observations
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
  
  ipmat <- createMAICInput(index = stackloss[c(1,10,21), ],
                           target = mtch.targets,
                           dictionary = mtch.dict,
                           matching.variables = c("airflow.mean",
                                                  "airflow.var",
                                                  "acidconc.prop",
                                                  "watertemp.median")
                           )
  expect_warning(wts <- maicWeight(ipmat, 
                                   control = list(reltol = 1e-12,
                                          maxit = 5000)),
                 "Optimisation has not converged")
})