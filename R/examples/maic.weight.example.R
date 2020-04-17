target <- c("Air.Flow" = 60,
            "Water.Temp" = 21,
            "Prop.Acid.Conc.LT.90" = 0.7,
            "min.air.flow" = 55)

stackloss$match.conc.lt.90 <- 
  ifelse(stackloss$Acid.Conc. < 90, 1, 0)

dict <- data.frame(
  "match.id" = 
    c("airflow", "watertemp", 
      "acidconc", "min.airflow"),
  "target.variable" = 
    c("Air.Flow", "Water.Temp", 
      "Prop.Acid.Conc.LT.90", "min.air.flow"),
  "index.variable" = 
    c("Air.Flow", "Water.Temp", 
      "match.conc.lt.90", "Air.Flow"),
  "match.type" = 
    c("mean", "mean", "proportion", "min"), 
  stringsAsFactors = FALSE)

weightObj <- maicMatching(
  index = stackloss,
  target = target,
  dictionary = dict,
  matching.variables = 
    c("airflow", "watertemp", 
      "acidconc", "min.airflow"))
