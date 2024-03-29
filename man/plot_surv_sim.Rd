% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation_plots.R
\name{plot_surv_sim}
\alias{plot_surv_sim}
\title{Plot Survival Function Simulations}
\usage{
plot_surv_sim(surv_sim, init_age, target_year, level = 95, years = NULL)
}
\arguments{
\item{surv_sim}{3D array of forecasted survival functions with survival time
(on the rows) and calendar year or cohort (on the columns) and simulation
number (3rd dimension)}

\item{init_age}{initial age for which `surv_sim` was calculated at}

\item{target_year}{year for which the survival function is plotted for}

\item{level}{desired confidence level with 95\% as default}

\item{years}{optional vector of years for `surv_sim`. If not supplied, then the column names
of `surv_sim` will be preserved}
}
\description{
Plots simulated survival functions for a desired cohort.
}
\examples{
# generate simulated rates with 'StMoMo'
# install and load 'StMoMo' if the package is not loaded
\dontrun{
# fitting lee carter model on ages 55:89
AUS_StMoMo <- StMoMoData(mortality_AUS_data, series = "male")
LC <- lc(link = "logit") # lee carter model
AUS_Male_Ini_Data <- central2initial(AUS_StMoMo)
ages_fit <- 55:89
wxy <- genWeightMat(ages = ages_fit, years = AUS_Male_Ini_Data$years, clip = 3)
LC_fit <- fit(LC, data = AUS_Male_Ini_Data, ages.fit = ages_fit, wxt = wxy)

# simulating rates for next 100 years
set.seed(1234)
n_sim <- 10
LC_sim <- simulate(LC_fit, nsim = n_sim, h = 100)

# using kannisto method to complete rates
young_ages <- LC_sim$ages # 55:89
old_ages <- 90:130
ages <- c(young_ages, old_ages)

kannisto_sim <- complete_old_age(rates = LC_sim$rates, ages = young_ages,
                                 old_ages = old_ages, control = list(fitted_ages = 80:89),
                                 method = "kannisto", type = "central")

# create period survival function for individual aged 55
surv_sim <- rate2survival(kannisto_sim, ages, from = "central")

plot_surv_sim(surv_sim, 55, 2050)
}
}
