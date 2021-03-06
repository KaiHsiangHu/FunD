% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FunD.R
\name{estimateAUC}
\alias{estimateAUC}
\title{Compute functional diversity by area under curve with particular sample coverages}
\usage{
estimateAUC(
  data,
  distM,
  datatype = "abundance",
  q = c(0, 1, 2),
  level = NULL,
  nboot = 50,
  conf = 0.95,
  tau = NULL
)
}
\arguments{
\item{data}{a matrix/data.frame of species abundances/incidences data.\cr}

\item{distM}{a pair wise distance matrix for all pairs of observed species in the pooled assemblage.\cr}

\item{datatype}{data type of input data: individual-based abundance data (\code{datatype = "abundance"}),
or species by sampling-units incidence frequencies (\code{datatype = "incidence_freq"}), default is \code{"abundance"}. \cr}

\item{q}{a sequence of nonnegative integers specifying the diversity orders of FD. Default is \code{c(0,1,2)}. \cr}

\item{level}{a positive sequence < 1 specifying a particular values of sample coverages.
If \code{NULL},then \code{level} will be chosen as the minimum coverage of all sites after extrapolating each site to its double sample size. Default is \code{NULL}.}

\item{nboot}{a positive integer specifying the number of bootstrap replications. Enter 0 to skip bootstrap;
in this case, the caculation of standard errors and confidence intervals will be skipped. Default is 50.}

\item{conf}{a positive number < 1 specifying the level of confidence interval, default is 0.95. \cr}

\item{tau}{a sequence between 0 and 1 specifying tau for integrating area under curve. If \code{NULL}, \code{tau} = (0, 0.01, 0.02,..., 0.99, 1). Default is \code{NULL}.}
}
\value{
a table including the sample size, sample coverage,
method (Interpolated or Extrapolated), and diversity estimates with each \code{q} for the user-specified sample coverages. \cr\cr
}
\description{
\code{estimateAUC}: computes functional diversity(FD) by area under curve thorough several thresholds with particular user-specified levels of sample coverages.
See Chao et al. (2019) for pertinent background and methods.
}
\examples{
\donttest{
# Type (1) abundance data (treat incidence frequencies as abundances to save computation time.)
data(FunDdata.abu)
data <- FunDdata.abu$data
dij <-  FunDdata.abu$dij
out <- estimateAUC(data = data, distM = dij, datatype = "abundance")
# Type (2) incidence frequency data 
data(FunDdata.inc)
data <- FunDdata.inc$data
dij <-  FunDdata.inc$dij
out <- estimateAUC(data = data, distM = dij, datatype = "incidence_freq")
}
}
\references{
Chao, A., Chiu C.-H. and Jost, L. (2010). functional diversity measures based on Hill numbers. Philosophical Transactions of the Royal Society B., 365, 3599-3609.\cr\cr
Chao, A., Chiu, C.-H., Hsieh, T. C., Davis, T., Nipperess, D., and Faith, D. (2015) Rarefaction and extrapolation of functional diversity. Methods in Ecology and Evolution, 6, 380-388.\cr\cr
Hsieh, T. C. and Chao, A. (2017). Rarefaction and extrapolation: making fair comparison of abundance-sensitive functional diversity among multiple assemblages. Systematic Biology 66, 100-111.
}
