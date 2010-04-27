\name{calcSE}
\alias{calcSE}
\title{ Calculate Standard Error from Replicates }
\description{
 This function is used to calculate the standard error of simulation results.  It is used by several of the functions in MCE, \cr 
 and is called by specifying type="SE".
}
\usage{
calcSE(data, index)
}

\arguments{
  \item{data}{ vector.  Results from R replicates of a simulation. }
  \item{index}{ vector.  Indexes of data that the expectation is taken over. }
}

\value{
the standard error of the replication results.
}

\author{ Sebastien Haneuse }


\examples{
test <- rnorm(10000, mean=1, sd=0.5)
calcSE(test, seq(from=1, to =10000, by=1))

}

\keyword{univar  }