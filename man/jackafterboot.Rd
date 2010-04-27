\name{jackafterboot}
\alias{jackafterboot}

\title{ Jackknife-After-Bootstrap Method of MCE estimation }
\description{
 Efron and Tibshirani have introduced a method for estimating Monte Carlo error 
for the estimates generated from the bootstrap.  This method uses
ideas similar to a jacknife method.
}
\usage{
jackafterboot(data, B = 100, type = "", truth)
}

\arguments{
  \item{data}{ vector. replication results }
  \item{B}{ numeric (integer).  number of bootstrap samples to be taken }
  \item{type}{ character.  type of simulation result including "mean", "SE", "RE", and "PB" for mean, standard error, \cr 
  relative efficiency, and percent bias. }
  \item{truth}{numeric.  if type="PB" the true value.  }
}
\details{
 This function takes original data and perfoms the bootstrap as well as estimates
 the error in the bootstrap estimates. For an original data set with N points 
and B bootstrap samples, the first estimate is calculated using only the replicates 
that do not contain the first point of the data.  Then the second estimate is from 
only the replicates that do not contain the second point, and so on until all N 
points have been removed.  An estimate is calculated from these replications without 
the nth point.  Then the standard deviation from the N estimates should 
estimate Monte Carlo error.   
}
\value{
Jackknife-After-Bootstrap estimate of MCE
}
\references{ Efron, B. and R. Tibshirani (1993).  An Introduction to the Bootstrap.  New York: Chapman & Hall. }
\author{ Elizabeth Koehler }
\note{  
Efron notes that this technique runs into trouble if one point appears in all B bootstrap samples, but 
claims this is rare when N >= 10 and B >= 20.
}

\examples{
data1<-rnorm(1000)
jackafterboot(data1,100, type="mean")

}

\keyword{univar }

