#' Apply Functions Over Each Matrix of 3D Array
#'
#' Returns a 3D array of matrices obtained by applying a function over each matrix
#' of a 3D array.
#'
#' @param X
#' 3D array
#' @param FUN
#' the function to be applied
#'
#' @return
#' 3D array of matrices obtained by applying a function
#' @keywords internal
arr_apply <- function(X, FUN) {

  X_list <- lapply(seq(dim(X)[3]), function(i) matrix(X[, , i], nrow = dim(X)[1], ncol = dim(X)[2]))
  sapply(X_list, FUN, simplify = "array")
}

#' Generate Default Death Probabilities
#'
#' Generates default death probabilities to be used for simulating paths if
#' necessary.
#'
#' @param init_age
#' the initial age of the path
#' @param sex
#' character string denoting the gender of individuals, "F" for female and "M" for male
#' @param closure_age
#' maximum life span
#'
#' @return
#' a vector of 1-year death probabilities corresponding to the specified
#' initial and closure age
#'
#' @keywords internal
generate_default_qx <- function(init_age, sex = "F", closure_age = 130) {

  # Flagging errors ---------------------------------------------------------

  # init_age
  if (init_age < 55 | init_age > closure_age) {
    stop("initial age must be between 55 and the maximum age")
  }

  if (init_age != floor(init_age)) {
    stop("initial age must be an integer")
  }

  # sex
  if (sex != "F" & sex != "M") {
    stop("sex must be 'F' or 'M'")
  }

  # closure_age
  if (closure_age < 90 | closure_age != floor(closure_age)) {
    stop("maximum age must be an integer greater than 89")
  }


  # Implementation ----------------------------------------------------------

  # Generating default death probabilities
  young_ages <- 55:89

  if (sex == "F") {
    AUS_StMoMo <- StMoMo::StMoMoData(StMoMoExt::mortality_AUS_data, series = "female")
    rates_hist <- StMoMoExt::mortality_AUS_data$rate$female[as.character(young_ages), ]
  } else {
    AUS_StMoMo <- StMoMo::StMoMoData(StMoMoExt::mortality_AUS_data, series = "male")
    rates_hist <- StMoMoExt::mortality_AUS_data$rate$male[as.character(young_ages), ]
  }

  # Using M7 model to Forecast Rates
  M7 <- StMoMo::m7()
  AUS_Ini_Data <- StMoMo::central2initial(AUS_StMoMo)
  ages_fit <- young_ages
  wxy <- StMoMo::genWeightMat(ages = ages_fit, years = AUS_Ini_Data$years, clip = 3)
  M7_fit <- StMoMo::fit(M7, data = AUS_Ini_Data, ages.fit = ages_fit, wxt = wxy)
  M7_for <- forecast::forecast(M7_fit, h = 100)

  # Mortality Rate Completion with Kannisto Method
  old_ages <- 90:closure_age
  ages <- c(young_ages, old_ages)
  kannisto_hist <- complete_old_age(rates = rates_hist, ages = young_ages,
                                    old_ages = old_ages, control = list(fitted_ages = 80:89),
                                    closure_age = closure_age,
                                    method = "kannisto", type = "central")
  kannisto_for <- complete_old_age(rates = M7_for$rates, ages = young_ages,
                                   old_ages = old_ages, control = list(fitted_ages = 80:89),
                                   closure_age = closure_age,
                                   method = "kannisto", type = "central")

  # Combine Historical and Forecasted Rates
  kannisto_55_period <- cbind(kannisto_hist, kannisto_for)
  qx_period <- rate2rate(kannisto_55_period, from = "central", to = "prob")
  # Take year 2022
  death_probs <- period2cohort(qx_period, ages = ages, init_age = init_age)[, "2022"]

  return(death_probs)
}
