\name{mceCLT}
\alias{mceCLT}

\title{ Central Limit Theorem Method of MCE estimation }
\description{
 For some simulation-based estimates, the MCE estimate can be based on the central limit theorem and can be calculated as the standard 
deviation of the R estimates from a simulation divided by sqrt(R).  This method only applies
to simulation-based estimates that are generated as a mean of replication results.}
\usage{
mceCLT(data, type = "", truth = NULL)
}

\arguments{
  \item{data}{ vector. replication reults }
  \item{type}{ character.  type of simulation result of interest including "mean", "SE", "PB" and "RE" for mean, standard error, percent bias and relative efficiency }
  \item{truth}{ numeric.  if type="PB" the true value. }

}
\details{

}
\value{
CLT-based estimate of MCE
}
\references{  }
\author{ Elizabeth Koehler and Sebastien Haneuse }
\note{ This methods only applies to simulation-based estimates that are generated from a mean of replication results.
}
\seealso{  }
\examples{
test <- matrix(rnorm(10000, mean=1))
mceCLT(test, type="PB", truth=1)

%function(data, type="", truth=NULL)
%{
%  ##
%  R <- length(data)
%
%  ##
%  value <- NA
%  if(type == "mean") value <- sd(data) / sqrt(R)
%  if(type == "PB") value <- sd((data - truth) / truth * 100) / sqrt(R)
%  return(value)
%  }
}
    \keyword{univar}

