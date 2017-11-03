data("stackloss")
target <- c("Air.Flow" = 60,
            "Water.Temp" = 21,
            "Proportion.Acid.Conc.Less.Than.90" = 0.7)
stackloss$match.conc.lt.90 <- ifelse(stackloss$Acid.Conc. < 90, 1, 0)

dict <- data.frame("match.id" = c("airflow",
                                  "watertemp",
                                  "acidconc"),
                   "target.variable" = c("Air.Flow",
                                         "Water.Temp",
                                         "Proportion.Acid.Conc.Less.Than.90"),
                   "index.variable" = c("Air.Flow",
                                        "Water.Temp",
                                        "match.conc.lt.90"),
                   "match.type" = c("mean",
                                    "mean",
                                    "proportion"),
                   stringsAsFactors = FALSE)
ipmat <- createMAICInputMatrix(index = stackloss,
                               target = target,
                               dictionary = dict,
                               matching.variables = c("airflow",
                                                      "watertemp",
                                                      "acidconc"))
wt <- maicWeight(ipmat)
rcv <- reportCovariates(stackloss,
                        target,
                        dict,
                        matching.variables = c("airflow",
                                               "watertemp",
                                               "acidconc"),
                        wt)
