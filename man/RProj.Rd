\name{RProj}
\alias{RProj}

\title{ The bootstrap-then-group implementation of the Bootstrap Grouping Prediction Plot for estimating R. }
\description{
     According to the CLT, the line relating MCE to 1/sqrt(R) will always pass through the origin.  
If a few of the points along this line can be estimated, a line of best fit with the intercept set 
to zero can be drawn.  After the line is estimated, an estimate of R can be generated for a given MCE.  
This function plots this line and reports an estimate of R.
}
\usage{
RProj(data, Rseq, B, bb = 1, type = "", truth = NULL, targetMCE, plot = FALSE, digits = 2, xAxis = NULL, yAxis = NULL, zAxis = NULL, ...)
}

\arguments{
 \item{data}{ vector.  replication results }
  \item{Rseq}{ vector.  group sizes, or R*, at which to estimate MCE in order to draw a prediction line}
  \item{B}{ numeric (integer).  when method="bootstrap", number of bootstrap replicates to use }
  \item{bb}{ numeric (integer).  numer of grouping permutations }
  \item{type}{ character.  type of simulation result of interest including "mean", "SE", "PB" and "RE" for mean, standard error, percent bias and relative efficiency }
  \item{truth}{ numeric.  if type="PB" the true value. }
  \item{targetMCE}{ vector or numeric.  MCE's for which the prediction of R is of interest }
  \item{plot}{ logical.  TRUE if plot is to be shown }
  \item{digits}{ numeric (integer).  number of digits to display for estimated MCE }
  \item{xAxis}{ vector.  sequence of points to be included on the xaxis.  Leaving as NULL the xaxis is calculated automatically.}
  \item{yAxis}{ vector.  sequence of points to be included on the yaxis.  Leaving as NULL the yaxis is calculated automatically.}
  \item{zAxis}{ vector.  sequence of points to be included on the top axis.  Leaving as NULL the zaxis is calculated automatically. }
  \item{\dots}{ other plotting parameters }
}
\details{
 
}
\value{
estimates of R from bootstrap grouping replication plot along with the plot, if specified

}
\references{  }
\author{ Sebastien Haneuse and Elizabeth Koehler }
\note{ 
This should be used cautiously when the Central Limit Theorem does not apply, or when the simulation-based estimate 
is not generated as an average of the replication results.  Note that the top axis presents the R that 
corresponds to the x-axis value of 1/sqrt(R) directly below.

}
\seealso{  \code{\link{mceCLT}}, \code{\link{mceBoot}}, \code{\link{RProj}}, \code{\link{groupbootMCE}} in package MCE }
\examples{
test <- as.matrix(cbind(rnorm(1000, mean=1, sd=0.5), rnorm(1000, mean=1, sd=1)))
RProj(test, Rseq=c(100, 500, 250), type="RE", B=1000, bb=5, plot=TRUE, targetMCE=c(0.25,1))

%function(data, Rseq, B, bb=1, type="", truth=NULL, targetMCE, plot=FALSE, digits=2, xAxis=NULL, yAxis=NULL, zAxis=NULL, ...)
%{
%  ##
%  R <- switch(type,
%              "mean"=length(data),
%              "PB"=length(data),
%              "SE"=length(data),
%              "RE"=nrow(data))

%  ##
%  ##
%  Ymce  <- matrix(NA, nrow=bb, ncol=length(Rseq))
%  Xmce  <- matrix(NA, nrow=bb, ncol=length(Rseq))
%  for(k in 1:bb){
%  bootRes <- matrix(NA, nrow=B, ncol=length(Rseq))
%  for(i in 1:length(Rseq))
%  {
%    subsets <- matrix(sample(1:R, Rseq[i]*B, replace=TRUE), ncol=B)
%    for(b in 1:B)
%    {
%      if(type == "mean") bootRes[b,i] <- calcEX(data, subsets[,b])
%      if(type == "PB")   bootRes[b,i] <- calcPB(data, subsets[,b], truth=truth)
%      if(type == "SE")   bootRes[b,i] <- calcSE(data, subsets[,b])
%      if(type == "RE")   bootRes[b,i] <- calcRE(data, subsets[,b])
%    }
%  }
%  ##
%  Ymce[k,]  <- apply(bootRes, 2, sd)
%  Xmce[k,]  <- 1 / sqrt(Rseq)
%  ##
%  }
%  beta1 <- lm(as.vector(Ymce) ~ as.vector(Xmce) - 1)$coef

%  ##
%  value <- round(as.numeric((targetMCE/beta1)^(-2)))
%  if(plot == TRUE)
%  {
%    ##
%    xRange <- c(0, ifelse(is.null(xAxis) == TRUE, 1/sqrt(min(Rseq)), max(xAxis)))
%    yRange <- c(0, ifelse(is.null(yAxis) == TRUE, max(Ymce), max(yAxis)))
%    plot(xRange, yRange, axes=FALSE, type="n", xlab=expression(1 / sqrt(R)), ylab="Monte Carlo error", ...)
%    ##
%    points(Xmce, Ymce, pch="X", cex=1.5)
%    if(is.null(xAxis) == TRUE) xAxis <- round(seq(from=0, to=xRange[2], length=5), digits=digits[1])
%    axis(1, at=xAxis)
%    if(is.null(yAxis) == TRUE) yAxis <- round(seq(from=0, to=yRange[2], length=5), digits=digits[1])
%    axis(2, at=yAxis)
%    if(is.null(zAxis) == TRUE) zAxis <- c(value, Rseq)
%    axis(3, at=1/sqrt(zAxis), labels=zAxis)
%    abline(0, beta1, col="red", lty=2)    
%    ##
%    segments(0, targetMCE,1/sqrt(value), targetMCE, col="blue", lwd=2)
%    segments(1/sqrt(value), 0, 1/sqrt(value), targetMCE, col="blue", lwd=2)
%    text(1/sqrt(value[1]) * 1.1, mean(yAxis[(length(yAxis)-c(0,1))])-seq(0,(length(value)-1)*.25*yAxis[2], by=.25*yAxis[2]), paste("R estimate =", round(sort(value, decreasing=FALSE), digits=digits[length(digits)]), sep=" "), pos=4)
%    box(bty="l")
%  }
%  return(value)
%  }
}
  \keyword{univar}
