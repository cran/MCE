\name{runsPlot}
\alias{runsPlot}

\title{ Runs Plot }
\description{
  Plots the simulation result as more replicates are added.
}
\usage{
runsPlot(data, type = "mean", thin = 10, burnin = 1, digits = 1, 
  truth = NULL, xAxis = NULL, yAxis = NULL, ...)
}

\arguments{
  \item{data}{ vector.  replication results }
  \item{type}{ character.  type of simulation result of interest including "mean", "SE", "PB" and "RE" for mean, standard error, \cr 
  percent bias and relative efficiency  }
  \item{thin}{ numeric. will plot every nth simulation result }
  \item{burnin}{ numeric.  number of replications to wait until plotting begins }
  \item{digits}{ numeric. number of digits used for y-axis labels. }
  \item{truth}{ numeric.  if type="PB" the true value. }
  \item{xAxis}{ vector.  sequence of points to be included on the xaxis.  Leaving as NULL the xaxis is calculated automatically. }
  \item{yAxis}{ vector.  sequence of points to be included on the yaxis.  Leaving as NULL the yaxis is calculated automatically. }
  \item{\dots}{ additional plot arguements }
}


\author{ Sebastien Haneuse and Elizabeth Koehler }

\examples{
test <- rnorm(10000, mean=1)
runsPlot(test, type="mean")

}

\keyword{univar }

