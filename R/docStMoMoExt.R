#' StMoMoExt.
#'
#' An extension of the 'StMoMo' R package to produce mortality rates
#' at older ages and risk-free survival probabilities. It also includes
#' functions to simulate life paths and output summary statistics.
#'
#' @name StMoMoExt
#' @docType package
#' @importFrom dplyr near
#' @importFrom fanplot fan
#' @importFrom forecast forecast
#' @importFrom graphics lines
#' @importFrom grDevices colorRampPalette rgb
#' @importFrom stats approxfun uniroot quantile pnorm
#' @importFrom utils head tail
NULL

#' Australian mortality data
#'
#' Mortality rates of Australians from the Human Mortality Database.
#'
#' \code{mortality_AUS_data} contains mortality rates of Australians
#' for the period 1970-2018 and for ages 0-110.
#'
#' @format An object of class \code{demogdata} with the following components:
#' \describe{
#'   \item{type}{type of object (mortality)}
#'   \item{age}{vector of ages}
#'   \item{year}{vector of years}
#'   \item{rate}{a list containing one or more rate matrices with age rows and year columns}
#'   \item{label}{label of the data.}
#'   \item{pop}{a list of the same form as \code{rate} but containing population numbers
#'   instead of demographic rates}
#' }
#' @source Human Mortality Database \url{http://www.mortality.org/}.
#' @references Human Mortality Database (2022). University of California,
#' Berkeley (USA), and Max Planck Institute for Demographic Research (Germany).
#' Available at \url{www.mortality.org}.
"mortality_AUS_data"
