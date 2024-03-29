% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rate_completion.R
\name{denuit_goderniaux}
\alias{denuit_goderniaux}
\title{Denuit and Goderniaux Method}
\usage{
denuit_goderniaux(
  rates,
  ages,
  old_ages,
  type = "prob",
  closure_age = 130,
  start_fit_age = 75,
  smoothing = FALSE,
  years = NULL
)
}
\arguments{
\item{rates}{vector, matrix or 3D array of mortality rates with age (on the rows) and
calendar year or cohort (on the columns) and simulation number (3rd dimension)}

\item{ages}{vector of ages for `rates`}

\item{old_ages}{vector of old ages for which `rates` is to be completed for}

\item{type}{character string representing the type of mortality rate for input and output.
Takes the following values: "central" for central death rates,
"prob" for 1-year death probabilities, "force" for force of mortality}

\item{closure_age}{maximum life span}

\item{start_fit_age}{model is fitted to ages starting from this age}

\item{smoothing}{logical value indicating if smoothing is to be applied}

\item{years}{optional vector of years for `rates`. If not supplied, then the column names
of `rates` will be preserved}
}
\value{
completed mortality rates for all ages and old ages in the same format as
`rates`
}
\description{
Implements the Denuit and Goderniaux method of mortality rate completion for
old ages.
}
\examples{
AUS_male_rates <- mortality_AUS_data$rate$male
ages <- mortality_AUS_data$age # 0:110
old_ages <- 91:130

# first convert to death probabilities
AUS_male_qx <- rate2rate(AUS_male_rates, from = "central", to = "prob")

completed_qx <- denuit_goderniaux(AUS_male_qx, ages, old_ages, type = "prob")

# fit on ages 80:110 instead
completed_qx_from_80 <- denuit_goderniaux(
AUS_male_qx, ages, old_ages, type = "prob", start_fit_age = 80)
}
