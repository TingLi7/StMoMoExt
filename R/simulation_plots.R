
#' Plot Curtate Future Lifetime Simulations
#'
#' Plots simulated expected curtate future lifetime.
#'
#' @param exp_cfl_rates
#' matrix of simulated expected curtate future lifetime with simulation number
#' (on the rows) and calendar year or cohort (on the columns).
#' Should be generated from \code{\link{exp_cfl}}
#' @param years
#' vector of years for \code{exp_cfl_rates}
#' @param level
#' desired confidence level with 95\% as default
#'
#' @export
#'
#' @examples
#' # generate simulated rates with 'StMoMo'
#' # install and load 'StMoMo' if the package is not loaded
#'
#' \dontrun{
#' # fitting lee carter model on ages 55:89
#' AUS_StMoMo <- StMoMoData(mortality_AUS_data, series = "male")
#' LC <- lc(link = "logit") # lee carter model
#' AUS_Male_Ini_Data <- central2initial(AUS_StMoMo)
#' ages_fit <- 55:89
#' wxy <- genWeightMat(ages = ages_fit, years = AUS_Male_Ini_Data$years, clip = 3)
#' LC_fit <- fit(LC, data = AUS_Male_Ini_Data, ages.fit = ages_fit, wxt = wxy)
#'
#' # simulating rates for next 100 years
#' set.seed(1234)
#' n_sim <- 10
#' LC_sim <- simulate(LC_fit, nsim = n_sim, h = 100)
#'
#' # using kannisto method to complete rates
#' young_ages <- LC_sim$ages # 55:89
#' old_ages <- 90:130
#' ages <- c(young_ages, old_ages)
#'
#' # extracting necessary info from historical data and simulations
#' rates_hist <- mortality_AUS_data$rate$male[as.character(young_ages), ]
#' years_hist <- as.numeric(colnames(rates_hist))
#' years_sim <- LC_sim$years
#' years <- c(years_hist, years_sim)
#'
#' kannisto_sim <- complete_old_age(
#' rates = LC_sim$rates, ages = young_ages, old_ages = old_ages,
#' fitted_ages = 80:89, method = "kannisto", type = "central")
#' kannisto_hist <- complete_old_age(
#' rates = rates_hist, ages = young_ages, old_ages = old_ages,
#' fitted_ages = 80:89, method = "kannisto", type = "central")
#'
#' ################# USAGE BEGINS HERE ################
#'
#' # combining historical and simulations
#' kannisto_55_period <- combine_hist_sim(rates_hist = kannisto_hist, rates_sim = kannisto_sim)
#'
#' # working with cohort starting from age 55
#' kannisto_55 <- period2cohort(period_rates = kannisto_55_period, ages = ages)
#' kannisto_55_q <- rate2rate(kannisto_55, from = "central", to = "prob")
#'
#' exp_cfl_kannisto <- exp_cfl(qx = kannisto_55_q, ages = ages)
#' # Expected curtate future lifetime can only be computed for earlier (complete) cohorts
#' exp_cfl_kannisto_clean <- exp_cfl_kannisto[, as.character(1970:2043)]
#' plot_exp_cfl(exp_cfl_rates = exp_cfl_kannisto_clean, years = 1970:2043)
#' }
plot_exp_cfl <- function(exp_cfl_rates, years, level = 95) {


  # Flagging Errors ---------------------------------------------------------


  # exp_cfl_rates
  if (!is.numeric(exp_cfl_rates) | !is.matrix(exp_cfl_rates)) {
    stop("exp_cfl_rates must be a numeric matrix")
  }

  # years
  if (length(years) != NCOL(exp_cfl_rates)) {
    stop("length of years must be equal to number of columns of exp_cfl_rates")
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

  # level
  if (level > 100 | level < 0) {
    stop("level must be between 0 and 100")
  }


  # Implementation ----------------------------------------------------------

  # Calculating mean, upper and lower values of confidence interval for simulations
  # of expected curtate future lifetime
  exp_cfl_mean <- apply(exp_cfl_rates, 2, mean)
  exp_cfl_lower <- apply(exp_cfl_rates, 2, stats::quantile, 1 / 2 - level / 200)
  exp_cfl_upper <- apply(exp_cfl_rates, 2, stats::quantile, 1 / 2 + level / 200)

  # Computing x and y limits of plot
  plot_ylim <- range(exp_cfl_mean, exp_cfl_lower, exp_cfl_upper, na.rm = TRUE)
  plot_xlim <- range(years)

  # Initial plot of mean
  plot(x = years,
       y = exp_cfl_mean,
       xlim = plot_xlim,
       ylim = plot_ylim,
       type = "l",
       xlab = "Years",
       ylab = "Expected Curtate Future Lifetime (Years)")

  # Preparing fanplot parameters
  fan_col <- grDevices::colorRampPalette(c("grey60", grDevices::rgb(1, 1, 1)))
  fan_n <- 1

  # Adding confidence intervals
  fanplot::fan(rbind(exp_cfl_lower, exp_cfl_upper),
               data.type = "values",
               start = years[1],
               probs = c(1 / 2 - level / 200, 1 / 2 + level / 200),
               fan.col = fan_col, n.fan = fan_n + 1, ln = NULL)

  # Overlaying mean on top of confidence intervals
  graphics::lines(x = years, y = exp_cfl_mean)

}

#' Plot Survival Function Simulations
#'
#' Plots simulated survival functions for a desired cohort.
#'
#' @param surv_sim
#' 3D array of forecasted survival functions with survival time
#' (on the rows) and calendar year or cohort (on the columns) and simulation
#' number (3rd dimension)
#' @param init_age
#' initial age for which `surv_sim` was calculated at
#' @param target_year
#' year for which the survival function is plotted for
#' @param level
#' desired confidence level with 95\% as default
#' @param years
#' optional vector of years for `surv_sim`. If not supplied, then the column names
#' of `surv_sim` will be preserved
#'
#' @export
#'
#' @examples
#'
plot_surv_sim <- function(surv_sim, init_age, target_year, level = 95, years = NULL) {

  # Flagging Errors ---------------------------------------------------------

  # surv_sim
  if (!(is.array(surv_sim) & length(dim(surv_sim)) == 3)) {
    stop("simulated survival functions must be a 3D array")
  }

  if (!is.numeric(surv_sim)) {
    stop("survival function must be numeric")
  }

  if (any(surv_sim < 0, na.rm = T)) {
    stop("survival function must be non-negative")
  }

  if (any(surv_sim > 1, na.rm = T)) {
    stop("survival function must be less than or equal to 1")
  }

  # init_age
  if (init_age != floor(init_age)) {
    stop("initial age must be an integer")
  }

  # target_year
  if (!is.element(target_year, as.numeric(colnames(surv_sim)))) {
    stop("target year must be included in years for simulations")
  }

  if (!is.null(years) & !is.element(target_year, years)) {
    stop("target year must be included in years for simulations")
  }

  # level
  if (level > 100 | level < 0) {
    stop("level must be between 0 and 100")
  }

  # years
  if (!is.null(years)) {
    if (length(years) != NCOL(surv_sim)) {
      stop("length of years must be equal to number of columns of the simulated survival functions")
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

  # Extracting required variables from inputs
  surv_time <- as.numeric(rownames(surv_sim))
  surv_sim_cohort <- t(surv_sim[, as.character(target_year), ])

  # Calculating mean, upper and lower values of confidence interval for simulations
  # of survival function
  surv_mean <- apply(surv_sim_cohort, 2, mean)
  surv_lower <- apply(surv_sim_cohort, 2, stats::quantile, 1 / 2 - level / 200)
  surv_upper <- apply(surv_sim_cohort, 2, stats::quantile, 1 / 2 + level / 200)

  # Computing x and y limits of plot
  plot_ylim <- range(surv_mean, surv_lower, surv_upper, na.rm = TRUE)
  plot_xlim <- c(surv_time[1], utils::tail(surv_time, 1))

  # Generating new plot
  plot(NULL,
       xlim = plot_xlim,
       ylim = plot_ylim,
       xlab = "Survival Time (Years)",
       ylab = "Survival Probability",
       main = paste("Age", init_age, "in Year", target_year))

  # Preparing fanplot parameters
  fan_col <- grDevices::colorRampPalette(c("grey60", grDevices::rgb(1, 1, 1)))
  fan_n <- 1

  # Adding confidence intervals
  fanplot::fan(rbind(surv_lower, surv_upper),
               data.type = "values",
               start = surv_time[1],
               probs = c(1 / 2 - level / 200, 1 / 2 + level / 200),
               fan.col = fan_col, n.fan = fan_n + 1, ln = NULL)

  # Plot mean
  graphics::lines(x = surv_time, y = surv_mean, type = 'l')

}
