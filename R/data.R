#' Simulated trial dataset for iron deficiency anaemia
#' 
#' A simulated dataset representing a comparison of investigational 
#' intervention "B" against control intervention "A". To be used in conjunction
#' with "haemoglobin_aggregate_data" for demonstrating MAIC
#' @format ## `haemoglobin`
#' A data frame with 300 rows and 5 columns:
#' \describe{
#'   \item{arm}{Intervention arm}
#'   \item{sex}{"Female" or "Male"}
#'   \item{age}{Age at baseline (years)}
#'   \item{Hb_bl}{Serum haemoglobin at baseline (g/dL)}
#'   \item{Hb_delta}{Difference in serum haemoglobin between baseline and week 12 (g/dL)}
#'   }
"haemoglobin"

#' Aggregate data from two simulated iron deficiency anaemia trials
#' 
#' Summary baseline characteristics and outcomes from two simulated trials - 
#' "AB" measuring intervention "B" versus "A" per the data in "haemoglobin"
#' and "AC" measuring intervention "C" versus "A". Counterfactual outcomes are
#' provided for both trials - in "AB", replacing "B" with "C", and in "AC" 
#' replacing "C" with "B"
#' @format  ## "haemoglobin_aggregate_data"
#' A data frame with 10 rows and 11 columns:
#' \describe{
#'   \item{N}{Number of patients}
#'   \item{age_mean}{mean age at baseline (years)}
#'   \item{age_sd}{standard deviation of age at baseline (years)}
#'   \item{proportion_female}{proportion of patients female}
#'   \item{Hb_bl_mean}{mean serum haemoglobin at baseline (g/dL)}
#'   \item{Hb_bl_sd}{standard deviation of serum haemoglobin at baseline (g/dL)}
#'   \item{Hb_delta_mean}{mean difference in serum haemoglobin between baseline and week 12 (g/dL)}
#'   \item{Hb_delta_sd}{standard deviation of serum haemoglobin between baseline and week 12 (g/dL)}
#'   \item{id}{population identifier}
#'   \item{difference_Hb_delta_mean}{for rows representing paired study arms, difference in mean difference in serum haemoglobin between baseline and week 12 (g/dL)}
#'   \item{difference_Hb_delta_sd}{for rows representing paired study arms, standard deviation of difference in mean difference in serum haemoglobin between baseline and week 12 (g/dL)}
#'   }
"haemoglobin_aggregate_data"