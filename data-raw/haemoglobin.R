## code to prepare `haemoglobin` dataset goes here
set.seed(20230524)

# Study sizes
N_ab <- 300
N_ac <- 500

# Baseline characteristics
pFemale_ab <- 0.6
pFemale_ac <- 0.8

muAge_ab <- 45
sigmaAge_ab <- 10
muAge_ac <- 50
sigmaAge_ac <- 8

muHb_male_ab <- 9
sigmaHb_male_ab <- 1
muHb_female_ab <- 8.5
sigmaHb_female_ab <- 0.7

muHb_male_ac <- 8.8
sigmaHb_male_ac <- 0.8
muHb_female_ac <- 8.3
sigmaHb_female_ac <-0.5

# Outcomes - delta Hb
dHb_A <- 1
dHb_B <- 2.5
dHb_C <- 2

# Outcome modifiers
dHb_A_Female_lt50 <- -0.2
dHb_B_Female_lt50 <- -0.7
dHb_C_Female_lt50 <- -0.2
dHb_A_Female_gte50 <- 0
dHb_B_Female_gte50 <- -0.3
dHb_C_Female_gte50 <- 0

# Linear gradient below anaemic threshold. Let's say 12 female, 13 male
dHb_A_Female_grad_an <- 0.01
dHb_A_Male_grad_an <- 0.01
dHb_B_Female_grad_an <- 0.03
dHb_B_Male_grad_an <- 0.03
dHb_C_Female_grad_an <- 0.01
dHb_C_Male_grad_an <- 0.01

generate_blc <- function(N,
                         arms,
                         pFemale, 
                         muAge, sigmaAge,
                         muHb_F, sigmaHb_F,
                         muHb_M, sigmaHb_M){
  s <- data.frame(
    arm = rep(arms, each = N %/% length(arms)),
    sex = ifelse(runif(N) < pFemale,
                 "Female",
                 "Male"),
    age = rnorm(N, muAge, sigmaAge)
  )
  s[["Hb_bl"]] <- ifelse(
    s$sex == "Female",
    rnorm(N, muHb_F, sigmaHb_F),
    rnorm(N, muHb_M, sigmaHb_M)
  )
  s[["age_cat1"]] <- ifelse(
    s$age >= 50,
    "gte50",
    "lt50"
  )
  s[["target_delta_Hb"]] <- ifelse(
    s$sex == "Female",
    12 - s[["Hb_bl"]],
    13 - s[["Hb_bl"]]
  )
  s
}

generate_outcomes <- function(x, noise = 0.1){
  x[["Hb_delta"]] <- 
    x[["noise"]] +
    dHb_A * (x[["arm"]] == "A") +
    dHb_B * (x[["arm"]] == "B") +
    dHb_C * (x[["arm"]] == "C") +
    dHb_A_Female_lt50 * (x[["arm"]] == "A" & x[["sex"]] == "Female" & x[["age_cat1"]] == "lt50") +
    dHb_B_Female_lt50 * (x[["arm"]] == "B" & x[["sex"]] == "Female" & x[["age_cat1"]] == "lt50") +
    dHb_C_Female_lt50 * (x[["arm"]] == "C" & x[["sex"]] == "Female" & x[["age_cat1"]] == "lt50") +
    dHb_A_Female_gte50 * (x[["arm"]] == "A" & x[["sex"]] == "Female" & x[["age_cat1"]] == "gte50") +
    dHb_B_Female_gte50 * (x[["arm"]] == "B" & x[["sex"]] == "Female" & x[["age_cat1"]] == "gte50") +
    dHb_C_Female_gte50 * (x[["arm"]] == "C" & x[["sex"]] == "Female" & x[["age_cat1"]] == "gte50") +
    dHb_A_Female_grad_an * x[["target_delta_Hb"]] * (x[["arm"]] == "A" & x[["sex"]] == "Female") +
    dHb_A_Male_grad_an * x[["target_delta_Hb"]] * (x[["arm"]] == "A" & x[["sex"]] == "Male") +
    dHb_B_Female_grad_an * x[["target_delta_Hb"]] * (x[["arm"]] == "B" & x[["sex"]] == "Female") +
    dHb_B_Male_grad_an * x[["target_delta_Hb"]] * (x[["arm"]] == "B" & x[["sex"]] == "Male") +
    dHb_C_Female_grad_an * x[["target_delta_Hb"]] * (x[["arm"]] == "C" & x[["sex"]] == "Female") +
    dHb_C_Male_grad_an * x[["target_delta_Hb"]] * (x[["arm"]] == "C" & x[["sex"]] == "Male")
  x
}

studyAB <- generate_blc(N_ab,
                        c("A", "B"),
                        pFemale_ab,
                        muAge_ab, sigmaAge_ab,
                        muHb_male_ab, sigmaHb_male_ab,
                        muHb_female_ab, sigmaHb_female_ab)
studyAB[["noise"]] <- rnorm(nrow(studyAB), 0, 0.1)

studyAB <- generate_outcomes(studyAB)

# For the sake of interest, let's get the counterfactual for if pop AB were randomised to AC
studyAB_AC <- studyAB
studyAB_AC[["arm"]] <- ifelse(studyAB_AC[["arm"]] == "A",
                              "A",
                              "C")
studyAB_AC <- generate_outcomes(studyAB_AC)

studyAC <- generate_blc(N_ac,
                        c("A", "C"),
                        pFemale_ac,
                        muAge_ac, sigmaAge_ac,
                        muHb_male_ac, sigmaHb_male_ac,
                        muHb_female_ac, sigmaHb_female_ac)
studyAC[["noise"]] <- rnorm(nrow(studyAC), 0, 0.1)

studyAC <- generate_outcomes(studyAC)
#Likewise, let's get the counterfactual for if pop AC were randomised to AB
studyAC_AB <- studyAC
studyAC_AB[["arm"]] <- ifelse(studyAC_AB[["arm"]] == "A",
                              "A",
                              "B")
studyAC_AB <- generate_outcomes(studyAC_AB)

# Round observed variables
studyAB[["age"]] <- floor(studyAB[["age"]])
studyAB_AC[["age"]] <- floor(studyAB_AC[["age"]])
studyAC[["age"]] <- floor(studyAC[["age"]])
studyAC_AB[["age"]] <- floor(studyAC_AB[["age"]])

studyAB[["Hb_bl"]] <- round(studyAB[["Hb_bl"]], 2)
studyAB_AC[["Hb_bl"]] <- round(studyAB_AC[["Hb_bl"]], 2)
studyAC[["Hb_bl"]] <- round(studyAC[["Hb_bl"]], 2)
studyAC_AB[["Hb_bl"]] <- round(studyAC_AB[["Hb_bl"]], 2)

studyAB[["Hb_delta"]] <- round(studyAB[["Hb_delta"]], 2)
studyAB_AC[["Hb_delta"]] <- round(studyAB_AC[["Hb_delta"]], 2)
studyAC[["Hb_delta"]] <- round(studyAC[["Hb_delta"]], 2)
studyAC_AB[["Hb_delta"]] <- round(studyAC_AB[["Hb_delta"]], 2)

get_aggregate_statistics <- function(x){
  data.frame(
    N = nrow(x),
    age_mean = mean(x$age),
    age_sd = sd(x$age),
    proportion_female = mean(x$sex == "Female"),
    Hb_bl_mean = mean(x$Hb_bl),
    Hb_bl_sd = sd(x$Hb_bl),
    Hb_delta_mean = mean(x$Hb_delta),
    Hb_delta_sd = sd(x$Hb_delta)
  )
}

get_one_stat <- function(x, name){
  s <- get_aggregate_statistics(x)
  s[["id"]] <- name
  s
}

haemoglobin_aggregate_data <- mapply(
  get_one_stat,
  list(
    studyAB[studyAB$arm == "A", ],
    studyAB[studyAB$arm == "B", ],
    studyAB,
    studyAB_AC[studyAB_AC$arm == "C", ],
    studyAB_AC,
    studyAC[studyAC$arm == "A", ],
    studyAC[studyAC$arm == "C", ],
    studyAC,
    studyAC_AB[studyAC_AB$arm == "B", ],
    studyAC_AB
  ),
  c(
    "Study AB arm A",
    "Study AB arm B",
    "Study AB",
    "Study AB counterfactual arm C",
    "Study AB counterfactual AC",
    "Study AC arm A",
    "Study AC arm B",
    "Study AC",
    "Study AC counterfactual arm B",
    "Study AC counterfactual AB"
  ),
  SIMPLIFY = FALSE
)

haemoglobin_aggregate_data <- do.call(rbind, haemoglobin_aggregate_data)

haemoglobin_aggregate_data[c(3,5,8,10), "Hb_delta_mean"] <- NA
haemoglobin_aggregate_data[c(3,5,8,10), "Hb_delta_sd"] <- NA
haemoglobin_aggregate_data[["difference_Hb_delta_mean"]] <- NA
haemoglobin_aggregate_data[["difference_Hb_delta_sd"]] <- NA

haemoglobin_aggregate_data[3, "difference_Hb_delta_mean"] <- 
  haemoglobin_aggregate_data[2, "Hb_delta_mean"] -
  haemoglobin_aggregate_data[1, "Hb_delta_mean"]
haemoglobin_aggregate_data[3, "difference_Hb_delta_sd"] <- 
  haemoglobin_aggregate_data[2, "Hb_delta_sd"] +
  haemoglobin_aggregate_data[1, "Hb_delta_sd"]

haemoglobin_aggregate_data[5, "difference_Hb_delta_mean"] <- 
  haemoglobin_aggregate_data[4, "Hb_delta_mean"] -
  haemoglobin_aggregate_data[1, "Hb_delta_mean"]
haemoglobin_aggregate_data[5, "difference_Hb_delta_sd"] <- 
  haemoglobin_aggregate_data[4, "Hb_delta_sd"] +
  haemoglobin_aggregate_data[1, "Hb_delta_sd"]

haemoglobin_aggregate_data[8, "difference_Hb_delta_mean"] <- 
  haemoglobin_aggregate_data[7, "Hb_delta_mean"] -
  haemoglobin_aggregate_data[6, "Hb_delta_mean"]
haemoglobin_aggregate_data[8, "difference_Hb_delta_sd"] <- 
  haemoglobin_aggregate_data[7, "Hb_delta_sd"] +
  haemoglobin_aggregate_data[6, "Hb_delta_sd"]

haemoglobin_aggregate_data[10, "difference_Hb_delta_mean"] <- 
  haemoglobin_aggregate_data[9, "Hb_delta_mean"] -
  haemoglobin_aggregate_data[6, "Hb_delta_mean"]
haemoglobin_aggregate_data[10, "difference_Hb_delta_sd"] <- 
  haemoglobin_aggregate_data[9, "Hb_delta_sd"] +
  haemoglobin_aggregate_data[6, "Hb_delta_sd"]

haemoglobin <- studyAB[, c("arm", "sex", "age", "Hb_bl", "Hb_delta")]

usethis::use_data(haemoglobin, overwrite = TRUE)
usethis::use_data(haemoglobin_aggregate_data, overwrite = TRUE)