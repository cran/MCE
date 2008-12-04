\name{calcRE}
\alias{calcRE}
\title{ Calculate the Relative Efficiency }
\description{
    This function is used to calculate the relative efficiency of simulation results.  It is used by several of the functions in MCE, and is called by specifying type="RE".

}
\usage{
calcRE(data, index)
}

\arguments{
  \item{data}{ vector.  Results from R replicates of a simulation. }
  \item{index}{ vector.  Indexes of data that the expectation is taken over. }
}

\value{
    the relative efficiency of the replication results.
}
\references{  }
\author{ Sebastien Haneuse }

\seealso{  \code{\link{calcEX}} }
\examples{
test <- as.matrix(cbind(rnorm(10000, mean=1, sd=0.5), rnorm(10000, mean=1, sd=1)))
calcRE(test, seq(from=1, to=10000, by=1) )

%function(data, index) sd(data[index,1]) / sd(data[index,2]) * 100
}

\keyword{univar  }
