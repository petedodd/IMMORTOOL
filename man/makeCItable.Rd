% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{makeCItable}
\alias{makeCItable}
\title{Make a table of results}
\usage{
makeCItable(list.of.results)
}
\arguments{
\item{list.of.results}{a named list of results from makeResultList, ie from ITBstats with N=#deaths appended}
}
\value{
a data frame of results for all elements in the list
}
\description{
For generating a results table from a list of individual study results
}
\details{
Generates a results table including CIs for analyses a)-d) of ITB stats as well as the fraction of deaths expected to occur in the control arm.
}
\examples{
set.seed(1234) #set random number seed

## create some results
ans.vanderVaart <- makeResultList(
                    mortality.times = c(7, 30, 90),
                    mortality.fracs = c(0.05, 0.21, 0.31),
                    treatment.times = c(6, 9, 13, 90),
                    treatment.fracs = c(0.09, 0.19, 0.28, 0.37),
                    N = 147, Tmax = 90, Tlandmark = 1,
                    plotfit = TRUE,
                    simulation.cohort.size = 5e3 #NOTE this is too small for convenient
                   )

## fake list duplicating above:
L <- list(ans.vanderVaart, ans.vanderVaart, ans.vanderVaart)
names(L) <- c('Study 1', 'Study 2', 'Study 3')

## create table from results list
tab <- makeCItable(L)
tab

}
\author{
Pete Dodd
}
