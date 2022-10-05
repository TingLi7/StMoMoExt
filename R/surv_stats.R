#' The Survival Function
#'
#' Distribution and quantile function  of the survival function.
#'
#' Linear interpolation is performed between the discrete points of the survival function
#'
#' @param surv_fun
#' vector of survival function
#' @param surv_time
#' vector of survival times
#' @param surv_prob
#' vector of survival probabilities
#'
#' @name surv
NULL

#' @rdname surv
#'
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
#' AUS_male_rates, ages, old_ages, method = "kannisto", type = "central",
#' control = list(fitted_ages = fitted_ages))
#'
#' all_ages <- 0:130
#' surv_func <- rate2survival(completed_rates, ages = all_ages, from = 'central', init_age = 55)
#'
#' # take vector of survival function (consider year 2017)
#' surv_func_2017 <- surv_func[, "2017"]
#'
#' # calculate probability of surviving 10 and 20 years
#' psurv(surv_func_2017, c(10, 20))
#'
#' # calculating the 80% and 95% quantile survival time
#' qsurv(surv_func_2017, c(0.8, 0.95))
psurv <- function(surv_fun, surv_time) {

  # Flagging Errors ---------------------------------------------------------

  # surv_fun
  if (!is.vector(surv_fun) | !is.numeric(surv_fun)) {
    stop("survival function must be a numeric vector")
  }

  if (any(surv_fun < 0, na.rm = T)) {
    stop("survival function must be non-negative")
  }

  if (any(surv_fun > 1, na.rm = T)) {
    stop("survival function must be less than or equal to 1")
  }

  # surv_time
  if (!is.vector(surv_time) | !is.numeric(surv_time)) {
    stop("survival time must be a numeric vector")
  }


  # Implementation ----------------------------------------------------------

  n <- length(surv_fun)

  # Checking bounds
  stopifnot(surv_time >= 0)

  # Linear interpolation
  surv_fun_approx <- stats::approxfun(0:(n - 1), surv_fun)

  result <- surv_fun_approx(surv_time)
  # Note: returns NA for surv_time >= n - 1

  return(result)

}

#' @rdname surv
#'
#' @export
qsurv <- function(surv_fun, surv_prob) {

  # Flagging Errors ---------------------------------------------------------

  # surv_fun
  if (!is.vector(surv_fun) | !is.numeric(surv_fun)) {
    stop("survival function must be a numeric vector")
  }

  if (any(surv_fun < 0, na.rm = T)) {
    stop("survival function must be non-negative")
  }

  if (any(surv_fun > 1, na.rm = T)) {
    stop("survival function must be less than or equal to 1")
  }

  # surv_prob

  if (!is.vector(surv_prob) | !is.numeric(surv_prob)) {
    stop("survival probabilities must be a numeric vector")
  }

  if (any(surv_prob < 0, na.rm = T)) {
    stop("survival probabilities must be non-negative")
  }

  if (any(surv_prob > 1, na.rm = T)) {
    stop("survival probabilities must be less than or equal to 1")
  }

  # surv_time


  # Implementation ----------------------------------------------------------


  n <- length(surv_fun)

  # Checking bounds
  stopifnot(surv_prob >= 0, surv_prob <= 1)

  # Creating list of linearly interpolated functions shifted by element of surv_prob
  surv_fun_list <- lapply(surv_prob, function(p) stats::approxfun(0:(n - 1), surv_fun - p))

  # Use uniroot to determine quantile
  sapply(surv_fun_list, function(fn) stats::uniroot(fn, c(0, (n - 1)))$root)

}
