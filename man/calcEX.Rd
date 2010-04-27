\name{calcEX}
\alias{calcEX}
\title{ Calculate Expectation from Replicates }
\description{
  This function is used to calculate the expectation of simulation results.  It is used by several of the functions in MCE, \cr 
  and is called by specifying type="mean".
}
\usage{
calcEX(data, index)
}
\arguments{
  \item{data}{ vector.  Results from R replicates of a simulation. }
  \item{index}{ vector.  Indexes of data that the expectation is taken over. }
}

\value{
the mean of the replication results.
}

\author{ Sebastien Haneuse}


\examples{
test <- as.matrix(cbind(rnorm(10000, mean=1, sd=0.5), rnorm(10000, mean=1, sd=1)))
calcRE(test, seq(1:10000))

}


\keyword{univar  }

