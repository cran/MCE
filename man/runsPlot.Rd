\name{runsPlot}
\alias{runsPlot}

\title{ Runs Plot }
\description{
  Plots the simulation result as a more replicates are added.
}
\usage{
runsPlot(data, type = "mean", thin = 10, burnin = 1, digits = 1, truth = NULL, xAxis = NULL, yAxis = NULL, ...)
}

\arguments{
  \item{data}{ vector.  replication results }
  \item{type}{ character.  type of simulation result of interest including "mean", "SE", "PB" and "RE" for mean, standard error, percent bias and relative efficiency  }
  \item{thin}{ numeric. will plot every nth simulation result }
  \item{burnin}{ numeric.  number of replications to wait until plotting begins }
  \item{digits}{  }
  \item{truth}{ numeric.  if type="PB" the true value. }
  \item{xAxis}{ vector.  sequence of points to be included on the xaxis.  Leaving as NULL the xaxis is calculated automatically. }
  \item{yAxis}{ vector.  sequence of points to be included on the yaxis.  Leaving as NULL the yaxis is calculated automatically. }
  \item{\dots}{ additional plot arguements }
}
\details{
  
}
\value{

}
\references{ }
\author{ Sebastien Haneuse and Elizabeth Koehler }
\note{  


}
\seealso{  }
\examples{
test <- matrix(rnorm(10000, mean=1))
runsPlot(test, type="mean")

%## The function is currently defined as
%function(data, type="mean", thin=10, burnin=1, digits=1, truth=NULL, xAxis=NULL, yAxis=NULL, ...)
%{
%  ##
%  ylab <- switch(type,
%                 mean="Expectation",
%                 PB="Percent Bias",
%                 SE="Standard Error",
%                 RE="Relative Efficiency")
%
%  ##
%  R <- switch(type,
%              mean=length(data),
%              PB=length(data),
%              SE=length(data),
%              RE=nrow(data))

%  ##
%  keep  <- seq(from=0, to=R, by=thin)[-1]
%  value <- rep(NA, length(keep))
%  for(i in 1:length(keep))
%  {
%    if(type == "mean") value[i] <- calcEX(data, c(1:keep[i]))
%    if(type == "PB") value[i] <- calcPB(data, c(1:keep[i]), truth=truth)
%    if(type == "SE") value[i] <- calcSE(data, c(1:keep[i]))
%    if(type == "RE") value[i] <- calcRE(data, c(1:keep[i]))
%  }

%  ##
%  xRange <- c(0, ifelse(is.null(xAxis) == TRUE, R, max(xAxis)))
%  yRange <- range(value[(keep > burnin)])
%  plot(xRange, yRange, axes=FALSE, type="n", xlab="Number of Replications, R", ylab=ylab, ...)
%  ##
%  if(is.null(xAxis) == TRUE) xAxis <- seq(from=0, to=R, length=5)
%  axis(1, at=xAxis)
%  if(is.null(yAxis) == TRUE) yAxis <- round(seq(from=yRange[1], to=yRange[2], length=5), digits=digits)
%  axis(2, at=yAxis)
%  box(bty="l")
%  lines(keep, value, lwd=2, col="red")
%  invisible()
%  }
}

\keyword{univar }

