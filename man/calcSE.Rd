\name{calcSE}
\alias{calcSE}
\title{ Calculate Standard Error from Replicates }
\description{
 This function is used to calculate the standard error of simulation results.  It is used by several of the functions in MCE, and is called by specifying type="SE".
}
\usage{
calcSE(data, index)
}

\arguments{
  \item{data}{ vector.  Results from R replicates of a simulation. }
  \item{index}{ vector.  Indexes of data that the expectation is taken over. }
}
\details{
}
\value{
the standard error of the replication results.
}
\references{  }
\author{ Sebastien Haneuse }

\seealso{  }
\examples{
test <- rnorm(10000, mean=1, sd=0.5)
calcSE(test, seq(from=1, to =10000, by=1))
%function(data, index) sd(data[index])
}

\keyword{univar  }