#' Convert death probabilities to survival function
#'
#' Converts mortality rates to the associated survival function.
#'
#' The survival function has survival time (starting from 0) on the rows.
#'
#' @param rates
#' vector, matrix or 3D array of mortality rates with age
#' (on the rows) and calendar year or cohort (on the columns) and simulation
#' number (3rd dimension)
#' @param ages
#' vector of ages for `rates`
#' @param from
#' character string representing the type of mortality rate to be converted
#' from. Takes the following values: "central" for central death rates, "prob"
#' for 1-year death probabilities, "force" for force of mortality.
#' @param init_age
#' initial age for which the survival function is to be calculated at. If not provided,
#' the survival function will be calculated for the smallest age supplied in `ages`
#' @param years
#' optional vector of years for `rates`. If not supplied, then the column names
#' of `rates` will be preserved
#'
#' @return
#' associated survival function as a 3D array if `rates` is an array, or as a matrix otherwise
#'
#' @export
#'
#' @examples
#' # consider the Kannisto completion method on male mortality rates
#' AUS_male_rates <- mortality_AUS_data$rate$male
#' ages <- mortality_AUS_data$age # 0:110
#' old_ages <- 91:130
#' fitted_ages <- 76:90
#'
#' completed_rates <- complete_old_age(
#' AUS_male_rates, ages, old_ages, method = "kannisto", type = "central",
#' fitted_ages = fitted_ages)
#'
#' # compute survival function of an individual aged 55
#' all_ages <- 0:130
#' surv_func <- rate2survival(completed_rates, ages = all_ages, from = 'central', init_age = 55)
rate2survival <- function(rates, ages, from = "prob", init_age = NULL, years = NULL) {

  # Flagging errors ---------------------------------------------------------

  # rates
  if (!is.vector(rates) & !is.matrix(rates) & !(is.array(rates) & length(dim(rates)) == 3)) {
    stop("rates must be a vector, 2D matrix or a 3D array")
  }

  if (!is.numeric(rates)) {
    stop("rates must be numeric")
  }

  if (any(rates < 0, na.rm = T)) {
    stop("rates must be non-negative")
  }

  if (from == "prob" & any(rates > 1, na.rm = T)) {
    stop("1-yr death probabilities must be less than or equal to 1")
  }

  # ages
  if (length(ages) != NROW(rates)) {
    stop("length of ages must be equal to number of rows of rates")
  }

  if (!is.vector(ages) | !all(ages == floor(ages))) {
    stop("ages must be a vector of integers")
  }

  if (is.unsorted(ages) | utils::tail(ages, 1) - ages[1] + 1 != length(ages)) {
    stop("ages must be increasing by 1 at each step")
  }

  if (any(ages < 0)) {
    stop("ages must be non-negative")
  }

  # from
  if (!is.element(from, c("central", "prob", "force"))) {
    stop("from must be 'central', 'prob' or 'force'")
  }



  # years
  if (!is.null(years)) {
    if (length(years) != NCOL(rates)) {
      stop("length of years must be equal to number of columns of rates")
    }

    if (!is.vector(years) | !all(years == floor(years))) {
      stop("years must be a vector of integers")
    }

    if (is.unsorted(years) | utils::tail(years, 1) - years[1] + 1 != length(years)) {
      stop("years must be increasing by 1 at each step")
    }

    if (any(years < 0)) {
      stop("years must be non-negative")
    }
  }

  # Implementation ----------------------------------------------------------

  # Converting to 1-year death probabilities
  qx <- rate2rate(rates, from, "prob")

  # Converting to 1-year survival probabilities

  if(is.null(init_age)) {
    init_age <- ages[1]
  } else if (!is.element(init_age, ages)) {
    stop("initial age must be in ages")
  }

  if(init_age == ages[1]) {
    px <- 1 - qx
  } else {
    px <- 1 - utils::tail(qx, ages[1] - init_age)
  }

  # Deal with R numerical error
  px[dplyr::near(px, 0)] <- 0

  # Calculating survival function
  if (is.vector(px)) {
    St <- rbind(1, matrix(cumprod(px)))
  } else if (is.matrix(px)) {
    St <- rbind(1, apply(px, 2, cumprod))
  } else if (is.array(px)) {
    St <- arr_apply(px, function(x) rbind(1, apply(x, 2, cumprod)))
  }

  rownames(St) <- as.character(0:(ages[length(ages)] - init_age + 1))
  colnames(St) <- if (is.null(years)) colnames(qx) else as.character(years)

  return(St)

}

#' Convert survival function to death probabilities
#'
#' Converts the survival function to the associated mortality rates.
#'
#' @param surv
#' vector, matrix or 3D array of the survival function with survival time starting from 0
#' (on the rows) and calendar year or cohort (on the columns) and simulation number (3rd dimension)
#' @param ages
#' vector of desired ages for the resulting 1-year death probabilities
#' @param to
#' character string representing the type of mortality rate to be converted
#' to. Takes the following values: "central" for central death rates, "prob"
#' for 1-year death probabilities, "force" for force of mortality.
#' @param years
#' optional vector of years for `surv`. If not supplied, then the column names
#' of `surv` will be preserved
#'
#' @return
#' associated mortality rates in the same format as `surv`
#' @export
#'
#' @examples
#' # create survival function for an individual aged 55
#' AUS_male_rates <- mortality_AUS_data$rate$male
#' ages <- mortality_AUS_data$age # 0:110
#' old_ages <- 91:130
#' fitted_ages <- 76:90
#'
#' completed_rates <- complete_old_age(
#' AUS_male_rates, ages, old_ages, method = "kannisto", type = "central", fitted_ages = fitted_ages)
#'
#' all_ages <- 0:130
#' surv_func <- rate2survival(completed_rates, ages = all_ages, from = 'central', init_age = 55)
#'
#' # note that we typically transfer from survival function to rate after
#' # converting a survival function from real world measure to risk-free measure
#'
#' # convert from P to Q measure survival function
#' # see the section on risk neutral probability
#' surv_func_Q <- survivalP2Q(surv_func, method = "wang", lambda = 1.5)
#'
#' # convert from survival function to mortality rates
#' central_rates_Q <- survival2rate(surv_func_Q, 55:130, to = 'central')
survival2rate <- function(surv, ages, to = "prob", years = NULL) {

  # Flagging Errors ---------------------------------------------------------
  # surv
  if (!is.vector(surv) & !is.matrix(surv) & !(is.array(surv) & length(dim(surv)) == 3)) {
    stop("survival function must be a vector, 2D matrix or a 3D array")
  }

  if (!is.numeric(surv)) {
    stop("survival function must be numeric")
  }

  if (any(surv < 0, na.rm = T)) {
    stop("survival function must be non-negative")
  }

  if (any(surv > 1, na.rm = T)) {
    stop("survival function must be less than or equal to 1")
  }

  # ages
  if (length(ages) != NROW(surv) - 1) {
    stop("length of ages and survival times do not match")
  }

  if (!is.vector(ages) | !all(ages == floor(ages))) {
    stop("ages must be a vector of integers")
  }

  if (is.unsorted(ages) | utils::tail(ages, 1) - ages[1] + 1 != length(ages)) {
    stop("ages must be increasing by 1 at each step")
  }

  if (any(ages < 0)) {
    stop("ages must be non-negative")
  }

  # to
  if (!is.element(to, c("central", "prob", "force"))) {
    stop("to must be 'central', 'prob' or 'force'")
  }

  # years
  if (!is.null(years)) {
    if (length(years) != NCOL(surv)) {
      stop("length of years must be equal to number of columns of the survival function")
    }

    if (!is.vector(years) | !all(years == floor(years))) {
      stop("years must be a vector of integers")
    }

    if (is.unsorted(years) | utils::tail(years, 1) - years[1] + 1 != length(years)) {
      stop("years must be increasing by 1 at each step")
    }

    if (any(years < 0)) {
      stop("years must be non-negative")
    }
  }

  # Implementation ----------------------------------------------------------

  px <- ifelse(utils::head(surv, -1) == 0, 0, utils::tail(surv, -1) / utils::head(surv, -1))
  qx <- 1 - px

  rates <- rate2rate(qx, from = "prob", to = to)


  if (!is.vector(rates)) {
    rownames(rates) <- as.character(ages)
    colnames(rates) <- if (is.null(years)) colnames(surv) else as.character(years)
  }

  return(rates)

}
