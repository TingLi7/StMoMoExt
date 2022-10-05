#' Combine Historical and Simulated Rates
#'
#' Combines a matrix of historical and 3D array of simulated rates.
#'
#' Connects the cohort/years of historical and simulated years together.
#' See \code{\link{plot_exp_cfl}} for usage
#'
#' @param rates_hist
#' Matrix of historical mortality rates with age (on the rows) and
#' calendar year or cohort (on the columns)
#' @param rates_sim
#' 3D array of simulated mortality rates with age (on the rows) and
#' calendar year or cohort (on the columns) and simulation number (3rd dimension)
#'
#' @return
#' 3D array of combined historical and simulated rates
#' @export
#'
#' @examples
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
#' rates_hist <- mortality_AUS_data$rate$male[as.character(young_ages), ]
#'
#' kannisto_sim <- complete_old_age(rates = LC_sim$rates, ages = young_ages,
#'                                  old_ages = old_ages, control = list(fitted_ages = 80:89),
#'                                  method = "kannisto", type = "central")
#' kannisto_hist <- complete_old_age(rates = rates_hist, ages = young_ages,
#'                                   old_ages = old_ages, control = list(fitted_ages = 80:89),
#'                                   method = "kannisto", type = "central")
#' ################# USAGE BEGINS HERE ################
#' # combining
#' kannisto_55_period <- combine_hist_sim(rates_hist = kannisto_hist,
#'                                        rates_sim = kannisto_sim)
#' }
combine_hist_sim <- function(rates_hist, rates_sim) {

  # Flagging errors ---------------------------------------------------------

  if (!is.matrix(rates_hist)) {
    stop("historical rates must be a 2D matrix")
  }

  if (!(is.array(rates_sim) & length(dim(rates_sim)) == 3)) {
    stop("simulated rates must be a 3D array")
  }

  if (dim(rates_hist)[1] != dim(rates_sim)[1]) {
    stop("historical and simulated rates must have an equal number of rows")
  }


  # Implementation ----------------------------------------------------------

  n_row = dim(rates_hist)[1]
  n_col = dim(rates_hist)[2] + dim(rates_sim)[2]
  n_sim = dim(rates_sim)[3]


  rates_all <- array(NA, dim = c(n_row, n_col, n_sim))
  for (i in 1:n_sim) {
    rates_all[,,i] = cbind(rates_hist, rates_sim[,,i])
  }

  rownames(rates_all) <- rownames(rates_hist)
  colnames(rates_all) <- c(colnames(rates_hist), colnames(rates_sim))

  rates_all
}
