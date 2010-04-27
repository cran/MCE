\name{calcRE}
\alias{calcRE}
\title{ Calculate the Relative Efficiency }
\description{
    This function is used to calculate the relative efficiency of simulation results.  It is used by several of the functions in \cr
     MCE, and is called by specifying type="RE".

}
\usage{
calcRE(data, index)
}

\arguments{
  \item{data}{ matrix.  Results from R replicates of a simulation, as two columns with efficiencies to compare. }
  \item{index}{ vector.  Indexes of data that the expectation is taken over. }
}

\value{
    the relative efficiency of the replication results.
}

\author{ Sebastien Haneuse }

\seealso{  \code{\link{calcEX}} }
\examples{
test <- as.matrix(cbind(rnorm(10000, mean=1, sd=0.5), rnorm(10000, mean=1, sd=1)))
calcRE(test, seq(from=1, to=10000, by=1) )


}

\keyword{univar  }
