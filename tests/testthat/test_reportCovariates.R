context("Reporting of covariates")
library(maic)

test_that("Presentation of weighted data", {
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
    "Water.Temp.median" = 22,
    "N" = 50
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
  sample.size.variable = c(
    "N", "N", "N", "N",
    "N", "N",
    "N", "N", "N"
  ),
  stringsAsFactors = FALSE
  )
  
  ipmat <- createMAICInput(index = stackloss,
                           target = mtch.targets,
                           dictionary = mtch.dict,
                           matching.variables = c("airflow.mean",
                                                  "airflow.sd",
                                                  "airflow.min",
                                                  "airflow.max",
                                                  "acidconc.prop",
                                                  "watertemp.median"))
  wts <- maicWeight(ipmat, control = list(reltol = 1e-12))
  
  expect_warning(rcv <- reportCovariates(stackloss,mtch.targets,
                   mtch.dict,
                   c("airflow.mean",
                     "airflow.sd",
                     "airflow.min",
                     "airflow.max",
                     "acidconc.prop",
                     "watertemp.median"),
                   wts),
                 "Chi-squared approximation may be incorrect")
  
  expect_equal(nrow(rcv), 6)
  expect_true(rcv["airflow.mean", "unadjusted.p.value"] < 0.05)
  expect_true(rcv["airflow.mean", "adjusted.p.value"] > 0.9)
  expect_true(rcv["airflow.sd", "unadjusted.p.value"] < 0.9)
  expect_true(rcv["airflow.sd", "adjusted.p.value"] > 0.9)
  expect_true(rcv["acidconc.prop", "unadjusted.p.value"] < 0.9)
  expect_true(rcv["acidconc.prop", "adjusted.p.value"] > 0.9)
  
  for (mtch in c("airflow.mean",
                 "airflow.sd",
                 "airflow.min",
                 "airflow.max",
                 "acidconc.prop",
                 "watertemp.median")){
    expect_true(abs(rcv[mtch, "adjusted.delta"]) <= abs(rcv[mtch, "unadjusted.delta"]))
  }
  
  expect_warning(rcv <- reportCovariates(stackloss,
                                         mtch.targets,
                                         mtch.dict,
                                         c("airflow.mean",
                                           "airflow.sd",
                                           "airflow.var",
                                           "airflow.min",
                                           "airflow.max",
                                           "acidconc.prop",
                                           "watertemp.10tile", 
                                           "watertemp.90tile",
                                           "watertemp.median"),
                                         wts),
                 "Chi-squared approximation may be incorrect")
  expect_equal(nrow(rcv), 9)
  
  ##### Subgrouping #####
  mtch.targets <- list("Air.Flow.mean" = 65,
                       "Proportion.Acid.Conc.lt.90" = 0)
  ipmat <- createMAICInput(index = stackloss,
                           target = mtch.targets,
                           dictionary = mtch.dict,
                           matching.variables = c("airflow.mean",
                                                  "acidconc.prop")
  )
  wts <- maicWeight(ipmat, control = list(reltol = 1e-12))
  # Expect no error
  expect_error(rcv <- reportCovariates(stackloss,mtch.targets,
                                       mtch.dict,
                                       c("acidconc.prop"),
                                       wts),
               NA)
  
  mtch.targets <- list("Air.Flow.mean" = 65,
                       "Proportion.Acid.Conc.lt.90" = 1)
  ipmat <- createMAICInput(index = stackloss,
                           target = mtch.targets,
                           dictionary = mtch.dict,
                           matching.variables = c("airflow.mean",
                                                  "acidconc.prop")
  )
  wts <- maicWeight(ipmat, control = list(reltol = 1e-12))
  # Expect no error
  expect_error(rcv <- reportCovariates(stackloss,mtch.targets,
                                       mtch.dict,
                                       c("acidconc.prop"),
                                       wts),
               NA)
})