% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/period_cohort.R
\name{period2cohort}
\alias{period2cohort}
\title{Convert from Period to Cohort Rates}
\usage{
period2cohort(period_rates, ages, init_age = NULL)
}
\arguments{
\item{period_rates}{matrix or 3D array of mortality rates with age (on the rows) and calendar
year (on the columns) and simulation number (3rd dimension)}

\item{ages}{vector of ages for `period_rates`}

\item{init_age}{initial age for which the cohort is to be considered. If not provided,
the cohort will correspond to the smallest age supplied in `ages`}
}
\value{
associated cohort mortality rates as a 3D array if `period_rates` is an array,
or as a matrix otherwise
}
\description{
Converts from period to cohort mortality rates
}
\examples{
period_rates <- mortality_AUS_data$rate$male
ages <- mortality_AUS_data$age # 0:110

# convert to rates for cohort aged 55
cohort_rates_55 <- period2cohort(period_rates, ages, init_age = 55)
}
