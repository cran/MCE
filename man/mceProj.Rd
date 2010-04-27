\name{mceProj}
\alias{mceProj}
\title{ The bootstrap-then-group implementation of the Bootstrap Grouping Prediction Plot for estimating MCE. }
\description{
    According to the CLT, the line relating MCE to 1/sqrt(R) will always pass through the origin.  
If a few of the points along this line can be estimated, a line of best fit with the intercept set 
to zero can be drawn.  After the line is estimated, an estimate of MCE can be generated for a given R.  
This function plots this line and reports an estimate of MCE.
}
\usage{
mceProj(data, Rseq, B, bb = 1, targetR = NULL, type = "", truth = 
  NULL, plot = FALSE, digits = 2, xAxis = NULL, yAxis = NULL, 
  zAxis = NULL, ...)
}
\arguments{
  \item{data}{ vector.  replication results }
  \item{Rseq}{ vector.  group sizes, or R*, at which to estimate MCE in order to draw a prediction line}
  \item{B}{ numeric (integer).  when method="bootstrap", number of bootstrap replicates to use }
  \item{bb}{ numeric (integer).  numer of grouping permutations }
  \item{targetR}{ vector or numeric (integer).  R's for which the prediction of MCE is of interest }
  \item{type}{ character.  type of simulation result of interest including "mean", "SE", "PB" and "RE" for mean, standard error, \cr
   percent bias and relative efficiency }
  \item{truth}{ numeric.  if type="PB" the true value. }
  \item{plot}{ logical.  TRUE if plot is to be shown }
  \item{digits}{ numeric (integer).  number of digits to display for estimated MCE }
  \item{xAxis}{ vector.  sequence of points to be included on the xaxis.  Leaving as NULL the xaxis is calculated automatically.}
  \item{yAxis}{ vector.  sequence of points to be included on the yaxis.  Leaving as NULL the yaxis is calculated automatically.}
  \item{zAxis}{ vector.  sequence of points to be included on the top axis.  Leaving as NULL the zaxis is calculated automatically. }
  \item{\dots}{ arguments passed to the plot function }
  
}

\value{
  estimates of MCE from bootstrap grouping replication plot along with the plot, if specified

}

\author{ Sebastien Haneuse and Elizabeth Koehler }
  \note{ 
 This should be used cautiously when the Central Limit Theorem does not apply, or when the simulation-based estimate 
is not generated as an average of the replication results.  Note that the top axis presents the R that 
corresponds to the x-axis value of 1/sqrt(R) directly below.
}
\seealso{ \code{\link{mceCLT}}, \code{\link{mceBoot}}, \code{\link{RProj}}, \code{\link{groupbootMCE}} in package MCE.}

\examples{
test <- as.matrix(cbind(rnorm(1000, mean=1, sd=0.5), rnorm(1000, mean=1, sd=1)))
mceProj(test, Rseq=c(100, 500, 250), type="RE", B=1000, bb=5, plot=TRUE,
  targetR=c(5000,2500))


}
 \keyword{univar  }
