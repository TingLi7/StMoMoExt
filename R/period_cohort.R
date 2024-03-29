#' Convert from Period to Cohort Rates
#'
#' Converts from period to cohort mortality rates
#'
#' @param period_rates
#' matrix or 3D array of mortality rates with age (on the rows) and calendar
#' year (on the columns) and simulation number (3rd dimension)
#' @param ages
#' vector of ages for `period_rates`
#' @param init_age
#' initial age for which the cohort is to be considered. If not provided,
#' the cohort will correspond to the smallest age supplied in `ages`
#'
#' @return
#' associated cohort mortality rates as a 3D array if `period_rates` is an array,
#' or as a matrix otherwise
#' @export
#'
#' @examples
#' period_rates <- mortality_AUS_data$rate$male
#' ages <- mortality_AUS_data$age # 0:110
#'
#' # convert to rates for cohort aged 55
#' cohort_rates_55 <- period2cohort(period_rates, ages, init_age = 55)
period2cohort <- function(period_rates, ages, init_age = NULL) {


  # Flagging errors ---------------------------------------------------------

  # period_rates

  if (!is.matrix(period_rates) & !(is.array(period_rates) & length(dim(period_rates)) == 3)) {
    stop("period rates must be a 2D matrix or a 3D array")
  }

  if (!is.numeric(period_rates)) {
    stop("period rates must be numeric")
  }

  if (any(period_rates < 0, na.rm = T)) {
    stop("period rates must be non-negative")
  }


  # ages
  if (length(ages) != NROW(period_rates)) {
    stop("length of ages must be equal to number of rows of period rates")
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


  # Implementation ----------------------------------------------------------

  if(is.null(init_age)) {
    init_age <- ages[1]
  } else if (!is.element(init_age, ages)) {
    stop("initial age must be in ages")
  }

  # Extract relevant rates starting from init_age
  if(init_age != ages[1]) {
    period_rates <- utils::tail(period_rates, ages[1] - init_age)
  }

  p2c_mat <- function(p_rates) {

    # Convert vector to matrix if necessary
    p_mat <- as.matrix(p_rates)

    n_row <- nrow(p_mat)
    n_col <- ncol(p_mat)

    c_rates <- matrix(NA, nrow = n_row, ncol = n_col)

    # Filling out cohort rates
    for (i in 1:min(n_row, n_col)) {
      c_rates[i, ] <- c(p_mat[i, i:n_col], rep(NA, i - 1))
    }

    c_rates
  }


  if (is.matrix(period_rates)) {
    cohort_rates <- p2c_mat(period_rates)
  } else if (is.array(period_rates)) {
    cohort_rates <- arr_apply(period_rates, p2c_mat)
  }

  dimnames(cohort_rates) <- dimnames(period_rates)

  cohort_rates

}


#' Convert from Cohort to Period Rates
#'
#' Converts from cohort to period mortality rates
#'
#' @param cohort_rates
#' matrix or 3D array of mortality rates with age (on the rows) and cohort
#' (on the columns) and simulation number (3rd dimension)
#'
#' @return
#' associated period mortality rates as a 3D array if `cohort_rates` is an array,
#' or as a matrix otherwise
#' @export
#'
#' @examples
#' AUS_male_rates <- mortality_AUS_data$rate$male
#' ages <- mortality_AUS_data$age # 0:110
#' old_ages <- 91:130
#' fitted_ages <- 76:90
#'
#' completed_rates <- complete_old_age(
#' AUS_male_rates, ages, old_ages, method = "kannisto", type = "central",
#' control = list(fitted_ages = fitted_ages))
#'
#' # suppose these are rates for cohort starting at age 60
#' cohort_rates_60 <- completed_rates[as.character(60:130), ]
#' period_rates <- cohort2period(cohort_rates_60)
cohort2period <- function(cohort_rates) {

  # Flagging errors ---------------------------------------------------------

  # cohort_rates

  if (!is.matrix(cohort_rates) & !(is.array(cohort_rates) & length(dim(cohort_rates)) == 3)) {
    stop("cohort rates must be a 2D matrix or a 3D array")
  }

  if (!is.numeric(cohort_rates)) {
    stop("cohort rates must be numeric")
  }

  if (any(cohort_rates < 0, na.rm = T)) {
    stop("cohort rates must be non-negative")
  }


  # Implementation ----------------------------------------------------------

  c2p_mat <- function(c_rates) {

    # Convert vector to matrix if necessary
    c_mat <- as.matrix(c_rates)

    n_row <- nrow(c_mat)
    n_col <- ncol(c_mat)

    p_rates <- matrix(NA, nrow = n_row, ncol = n_col)

    # Filling out period rates
    for (i in 1:min(n_row, n_col)) {
      p_rates[i, ] <- c(rep(NA, i - 1), c_mat[i, 1:(n_col + 1 - i)])
    }

    p_rates

  }

  if (is.vector(cohort_rates) | is.matrix(cohort_rates)) {
    period_rates <- c2p_mat(cohort_rates)
  } else if (is.array(cohort_rates)) {
    period_rates <- arr_apply(cohort_rates, c2p_mat)
  }

  dimnames(period_rates) <- dimnames(cohort_rates)

  period_rates

}
