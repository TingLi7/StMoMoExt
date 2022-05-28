#' Summarise Curtate Future Lifetime Statistics
#'
#' Produces expected curtate future lifetime for a 3D array of 1-year
#' death probabilities.
#'
#' @param qx
#' vector, matrix or 3D array of 1-year death probabilities with age
#' (on the rows) and calendar year or cohort (on the columns) and simulation
#' number (3rd dimension)
#' @param ages
#' vector of ages for `qx`
#' @param init_age
#' initial age for which the curtate future lifetime is to be calculated at. If not provided,
#' the summary statistics will be calculated for the smallest age supplied in `ages`
#' @param years
#' optional vector of years for `qx`. If not supplied, then the column names
#' of `qx` will be preserved
#'
#' @return
#' expected curtate future lifetime as a matrix if `qx` is a 3D array
#' with simulation number (on the rows) and calendar year or cohort (on the columns).
#' Returns a vector otherwise
#'
#' @export
#'
#' @examples
#' # complete rates using kannisto method
#' AUS_male_rates <- mortality_AUS_data$rate$male
#' ages <- mortality_AUS_data$age # 0:110
#' old_ages <- 91:130
#' fitted_ages <- 76:90
#'
#' completed_rates <- complete_old_age(
#' AUS_male_rates, ages, old_ages, method = "kannisto", type = "central", fitted_ages = fitted_ages)
#'
#' # convert to qx
#' completed_qx <- rate2rate(completed_rates, from = "central", to = "prob")
#' all_ages <- 0:130
#'
#' # expected curtate future lifetime using period rates for individual aged 0 and 60
#' ecfl_0 <- exp_cfl(completed_qx, all_ages)
#' ecfl_60 <- exp_cfl(completed_qx, all_ages, init_age = 60)
exp_cfl <- function(qx, ages, init_age = NULL, years = NULL) {

  # Flagging Errors ---------------------------------------------------------

  # qx
  if (!is.vector(qx) & !is.matrix(qx) & !(is.array(qx) & length(dim(qx)) == 3)) {
    stop("qx must be a vector, 2D matrix or a 3D array")
  }

  if (!is.numeric(qx)) {
    stop("qx must be numeric")
  }

  if (any(qx < 0, na.rm = T)) {
    stop("qx must be non-negative")
  }

  if (any(qx > 1, na.rm = T)) {
    stop("qx must be less than or equal to 1")
  }

  # ages
  if (length(ages) != NROW(qx)) {
    stop("length of ages must be equal to number of rows of qx")
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

  # years
  if (!is.null(years)) {
    if (length(years) != NCOL(qx)) {
      stop("length of years must be equal to number of columns of qx")
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

  # Converting to 1-year survival probabilities
  if(is.null(init_age)) {
    init_age <- ages[1]
  } else if (!is.element(init_age, ages)) {
    stop("invalid initial age")
  }

  if(init_age == ages[1]) {
    px <- 1 - qx
  } else {
    px <- 1 - utils::tail(qx, ages[1] - init_age)
  }

  # Calculating kpx
  if (is.vector(px)) {
    kpx <- matrix(cumprod(px))
  } else if (is.matrix(px)) {
    kpx <- apply(px, 2, cumprod)
  } else if (is.array(px)) {
    kpx <- arr_apply(px, function(x) apply(x, 2, cumprod))
  }

  # Changing dim names
  stopifnot(dim(px) == dim(kpx))
  colnames(kpx) <- colnames(px)
  k <- 1:nrow(kpx)
  rownames(kpx) <- as.character(k)

  exp_cfl_mat <- function(kpx_mat) {
    # Expected curtate future lifetime
    # Assumes kpx has been given until terminal age

    return(apply(kpx_mat, 2, sum))
  }

  # kpx should be matrix or array, note that is.array(A) = TRUE where A is matrix
  stopifnot(is.array(kpx))
  if (is.matrix(kpx)) {
    result <- exp_cfl_mat(kpx)

    return(result)
  } else {
    result <- arr_apply(kpx, exp_cfl_mat)
    rownames(result) <- if (is.null(years)) colnames(qx) else as.character(years)

    return(t(result))
  }

}
