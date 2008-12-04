\name{mceBoot}
\alias{mceBoot}
\title{ Bootstrap-After-Bootstrap estimate of MCE }
\description{
 Just like using the bootstrap to mimic multiple samples, the bootstrap 
can be used to mimic generating multiple sets of stochastic outcomes to
 estimate the MCE.
}
\usage{
mceBoot(data, B, type = "", truth = NULL)
}

\arguments{
  \item{data}{vector. replication reults }
  \item{B}{ numeric (integer). number of bootstrap replications to use. }
  \item{type}{ character.  type of simulation result of interest including "mean", "SE", "PB" and "RE" for mean, standard error, percent bias and relative efficiency }
  \item{truth}{ numeric.  if type="PB" the true value. }
}
\details{
  \item{1.}{Draw a sample with replacement of length R from the replication results.}
  \item{2.}{Calculate and record the bootstrap result as a function of the sampled replication results
 to generate a simulation-based estimate.}
  \item{3.}{Repeat steps 1 and 2, B times.}  
  \item{4.}{A standard deviation of these B simulation-based estimates can be used to obtain an estimate 
	of the Monte Carlo error associated with R replications.}  
}
\value{
 bootstrap-after-boostrap estimate of MCE
}
\references{  }
\author{ Elizabeth Koehler and Sebastien Haneuse }
\note{ Unlike CLT methods, this method is not limited to simulation-based estimates calculated as a mean
}
\seealso{ \code{\link{mceCLT}} }
\examples{
test <- matrix(rnorm(10000, mean=1))
mceBoot(test, type="PB", B=500, truth=1)

%function(data, B, type="", truth=NULL)
%{
%  ##
%  R <- switch(type,
%              "mean"=length(data),
%              "PB"=length(data),
%              "SE"=length(data),
%              "RE"=nrow(data))
%
%  ##
%  subsets <- matrix(sample(1:R, R*B, replace=TRUE), ncol=B)
%  bootRes <- rep(NA, B)
%  for(b in 1:B)
%  {
%    if(type == "mean") bootRes[b] <- calcEX(data, subsets[,b])
%    if(type == "PB")   bootRes[b] <- calcPB(data, subsets[,b], truth=truth)
%    if(type == "SE")   bootRes[b] <- calcSE(data, subsets[,b])
%    if(type == "RE")   bootRes[b] <- calcRE(data, subsets[,b])
%  }
%  value <- sd(bootRes)
%  return(value)
%  }
}
\keyword{univar  }

