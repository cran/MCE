\name{groupbootMCE}
\alias{groupbootMCE}
\title{ The group-then-bootstrap implementation of the Bootstrap Grouping Prediction Plot for estimating MCE}
\description{
    According to the CLT, the line relating MCE to 1/sqrt(R) will always pass through the origin.  
If a few of the points along this line can be estimated, a line of best fit with the intercept set 
to zero can be drawn.  After the line is estimated, an estimate of MCE can be generated for a given R.  
This function plots this line and reports an estimate of MCE.
}
\usage{
groupbootMCE(data, Rseq = c(100, 500), B = 20, bb = 1, targetR = NA, 
 method = "bootstrap", truth = NULL, plot = FALSE, 
 digits = 2, xAxis = NULL, yAxis = NULL, zAxis = NULL, type = "", ...)
}

\arguments{
  \item{data}{ vector.  replication results }
  \item{Rseq}{ vector.  group sizes, or R*, at which to estimate MCE in order to draw a prediction line}
  \item{B}{ numeric (integer).  when method="bootstrap", number of bootstrap replicates to use }
  \item{bb}{ numeric (integer).  numer of grouping permutations }
  \item{targetR}{ vector or numeric (integer).  R's for which the prediction of MCE is of interest }
  \item{method}{ character.  either "clt" or "bootstrap" indicating the method used to estimate MCE at R* }
  \item{plot}{ logical.  TRUE if plot is to be shown }
  \item{digits}{ numeric (integer).  number of digits to display for estimated MCE }
  \item{xAxis}{ vector.  sequence of points to be included on the xaxis.  Leaving as NULL the xaxis is calculated automatically.}
  \item{yAxis}{ vector.  sequence of points to be included on the yaxis.  Leaving as NULL the yaxis is calculated automatically.}
  \item{zAxis}{ vector.  sequence of points to be included on the top axis.  Leaving as NULL the zaxis is calculated automatically. }
  \item{type}{ character.  type of simulation result of interest including "mean", "SE", "PB" and "RE" for mean, standard error, \cr
  percent bias and relative efficiency }
  \item{truth}{ numeric.  if type="PB" the true value. }
  \item{\dots}{ other plotting parameters }
}
\details{
  The following steps are used to generate the points that are used to fit the line.
  \enumerate{
  \item To start generating points to use when plotting the prediction line, MCE estimates for simulations 
 with fewer than R replications are used.  In order to recover the distribution that would have been 
observed if a simulation with R* were run many times, the R replications can be randomly broken into 
 groups of R* replicates.  

  \item Each of these groups will generate an estimate of MCE using the central limit theorem method, or the \cr
  bootstrap-after-bootstrap method that is then averaged to give an estimate of MCE when R is R*.
  
   \item The average of estimated MCE is plotted as a point with the x-coordinate 1/sqrt(R*), and 
y-coordinate as the estimate from the previous step. 
   \item The estimate from step 3 may have been affected by the way groups are assigned, 
so steps 1-3 are repeated bb times so that the estimate of MCE will be unaffected by group assignments.  
In other words, the group assignments are permuted bb times.  
   \item Then, these points are used to plot a line through the origin by setting the intercept term to zero.  
   \item It may be desirable to use more than the estimates from a single R* to generate this line, and it is \cr
   recommended that multiple R* be used, or Rseq is a vector.
   \item After the CLT has been used to generate the line, an estimate of MCE for a specific R is possible.}
}

\value{
estimates of MCE from bootstrap grouping replication plot along with the plot, if specified
}

\author{ Elizabeth Koehler and Sebastien Haneuse}
\note{ 
 This should be used cautiously when the Central Limit Theorem does not apply, or when the simulation-based estimate 
is not generated as an average of the replication results.  Note that the top axis presents the R that 
corresponds to the x-axis value of 1/sqrt(R) directly below.
}
\seealso{ \code{\link{mceCLT}}, \code{\link{mceBoot}}, \code{\link{RProj}}, \code{\link{groupbootMCE}} in package MCE.}
\examples{
percentbias<-rnorm(10000,15, 4)
groupbootMCE(percentbias, Rseq=c(100,500,2000), B=20, bb=10, 
  targetR=c(15000, 10000), method="bootstrap", type="PB", 
  truth=15, plot=TRUE)
}

\keyword{univar  }

