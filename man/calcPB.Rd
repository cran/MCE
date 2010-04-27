\name{calcPB}
\alias{calcPB}
\title{ Calculate Percent Bias from Replicates }
\description{
   This function is used to calculate the percent bias of simulation results.  It is used by several of the functions in MCE, \cr 
   and is called by specifying type="PB".

}
\usage{
calcPB(data, index, truth)
}

\arguments{
  \item{data}{ vector.  Results from R replicates of a simulation. }
  \item{index}{ vector.  Indexes of data that the percent bias is taken over. }
  \item{truth}{ numeric.  The true value. }
}

\value{
 The percent bias estimated from the replicates.
}

\author{ Sebastien Haneuse}


\examples{
test <- rnorm(10000, mean=1)

calcPB(test, seq(from=1, to=10000, by=1), truth=1)

}
 \keyword{univar  }

