#' Survival Function Transformation
#'
#' Transforms the survival function from the real world P-measure to the risk-neutral
#' Q-measure according to the specified risk-neutral principles.
#'
#' The risk-neutral principles and their corresponding strings and valid lambda range
#' are as follows:
#' * "wang", \eqn{\lambda \ge 0}: Wang Transform,
#' * "ph", \eqn{\lambda \ge 1}: Proportional Hazard Transform
#' * "dp", \eqn{\lambda \ge 1}: Dual-power Transform
#' * "gp", \eqn{0 \le \lambda \le 1}: Gini Principle
#' * "dadp", \eqn{0 \le \lambda \le 1}: Denneberg's Absolute Deviation Principle
#' * "exp", \eqn{\lambda > 0}: Exponential Transform
#' * "log", \eqn{\lambda > 0}: Logarithmic Transform
#' * "canon", \eqn{\lambda > 0}: Univariate Canonical Valuation
#' * "esscher", \eqn{\lambda > 0}: Esscher Transform
#'
#' The first seven principles are distortion risk measures and act on the survival function.
#' The last two principles act on the probability density function.
#'
#' @param StP
#' vector, matrix or 3D array of the survival function under the P-measure with
#' survival time (on the rows) and calendar year or cohort (on the columns) and simulation
#' number (3rd dimension)
#' @param method
#' character string representing the distortion risk measure to be used. See "Details".
#' @param lambda
#' parameter associated with the distortion risk measure
#'
#' @return
#' the transformed survival function under the specified Q-measure
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
#' # convert from P to Q measure survival function
#' surv_func_Q <- survivalP2Q(surv_func, method = "wang", lambda = 1.5)
survivalP2Q <- function(StP, method, lambda) {

  # Flagging Errors ---------------------------------------------------------

  # StP
  if (!is.vector(StP) & !is.matrix(StP) & !(is.array(StP) & length(dim(StP)) == 3)) {
    stop("survival function must be a vector, 2D matrix or a 3D array")
  }

  if (!is.numeric(StP)) {
    stop("survival function must be numeric")
  }

  if (any(StP < 0, na.rm = T)) {
    stop("survival function must be non-negative")
  }

  if (any(StP > 1, na.rm = T)) {
    stop("survival function must be less than or equal to 1")
  }

  # method
  valid_methods <- c("wang", "ph", "dp", "gp", "dadp", "exp", "log", "canon", "esscher")
  if (!is.element(method, valid_methods)) {
    stop("invalid risk distortion method")
  }

  if (!is.numeric(lambda)) {
    stop("lambda must be numeric")
  }

  # Implementation -----------------------------------------------------------

  # Defining survival function distortion functions
  wang <- function(x, lam) {
    if (lam < 0) stop("invalid lambda value")
    return(1 - stats::pnorm(stats::qnorm(1 - x) - lam))
  }

  ph <- function(x, lam) {
    if (lam < 1) stop("invalid lambda value")
    return(x^(1 / lam))
  }

  dp <- function(x, lam) {
    if (lam < 1) stop("invalid lambda value")
    return(1 - (1 - x)^lam)
  }

  gp <- function(x, lam) {
    if (lam < 0 | lam > 1) stop("invalid lambda value")
    return((1 + lam) * x - lam * x^2)
  }

  dadp <- function(x, lam) {
    if (lam < 0 | lam > 1) stop("invalid lambda value")
    return(ifelse(x < 0.5, (1 + lam) * x, lam + (1 - lam) * x))
  }

  exp_tfm <- function(x, lam) {
    if (lam <= 0) stop("invalid lambda value")
    return((1 - exp(-lam * x)) / (1 - exp(-lam)))
  }

  log_tfm <- function(x, lam) {
    if (lam <= 0) stop("invalid lambda value")
    return(log(1 + lam * x) / log(1 + lam))
  }

  survival_type <- c("wang", "ph", "dp", "gp", "dadp", "exp", "log")
  pdf_type <- c("canon", "esscher")

  if (is.element(method, survival_type)) {
    if (method == "wang") distort <- wang
    else if (method == "ph") distort <- ph
    else if (method == "dp") distort <- dp
    else if (method == "gp") distort <- gp
    else if (method == "dadp") distort <- dadp
    else if (method == "exp") distort <- exp_tfm
    else if (method == "log") distort <- log_tfm

    return(distort(StP, lambda))
  } else if (is.element(method, pdf_type)) {
    # Risk-adjusted pdf is identical for univariate canonical valuation and
    # esscher transform

    if (lambda <= 0) stop("invalid lambda value")

    pdfP2Q <- function(StP_mat) {

      time <- 0:(nrow(StP_mat) - 1)
      # Calculating pdf
      ftP <- rbind(0, diff(1 - StP_mat))
      rownames(ftP) <- as.character(time)
      # ftP will sum up to 1 if survival function starts from 1 and ends at 0
      stopifnot(nrow(StP_mat) == nrow(ftP))

      # Canonical valuation

      num <- exp(lambda * time) * ftP
      denom <- apply(num, 2, sum)
      ftQ <- num / outer(rep(1, nrow(StP_mat)), denom)
      StQ_mat <- 1 - apply(ftQ, 2, cumsum)

      return(StQ_mat)
    }

    # Treating vectors, matrices and arrays differently
    if (is.vector(StP)) {
      StQ <- pdfP2Q(as.matrix(StP))
    } else if (is.matrix(StP)) {
      StQ <- pdfP2Q(StP)
    } else if (is.array(StP)) {
      StQ <- arr_apply(StP, pdfP2Q)
    }

    stopifnot(dim(StP) == dim(StQ))
    dimnames(StQ) <- dimnames(StP)

    return(StQ)

  } else {
    stop("invalid method type")
  }

}
