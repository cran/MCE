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
groupbootMCE(data, Rseq = c(100, 500), B = 20, bb = 1, targetR = NA, method = "bootstrap", truth = NULL, plot = FALSE, digits = 2, xAxis = NULL, yAxis = NULL, zAxis = NULL, type = "", ...)
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
  \item{type}{ character.  type of simulation result of interest including "mean", "SE", "PB" and "RE" for mean, standard error, percent bias and relative efficiency }
  \item{truth}{ numeric.  if type="PB" the true value. }
  \item{\dots}{ other plotting parameters }
}
\details{
  The following steps are used to generate the points that are used to fit the line.
  \item{1.}{To start generating points to use when plotting the prediction line, MCE estimates for simulations 
 with fewer than R replications are used.  In order to recover the distribution that would have been 
observed if a simulation with R* were run many times, the R replications can be randomly broken into 
 groups of R* replicates. } 

  \item{2.}{Each of these groups will generate an estimate of MCE using the central limit theorem method, or the bootstrap-after-bootstrap method 
that is then averaged to give an estimate of MCE when R is R*.}
   \item{3.}{The average of estimated MCE is plotted as a point with the x-coordinate 1/sqrt(R*), and 
y-coordinate as the estimate from the previous step. } 
   \item{4.}{The estimate from step 3 may have been affected by the way groups are assigned, 
so steps 1-3 are repeated bb times so that the estimate of MCE will be unaffected by group assignments.  
In other words, the group assignments are permuted bb times.}  
   \item{5.}{Then, these points are used to plot a line through the origin by setting the intercept term to zero.}  
   \item{6.}{It may be desirable to use more than the estimates from a single R* to generate this line, and it is recommended that multiple R* be used, or Rseq is a vector.} 
   \item{7.}{After the CLT has been used to generate the line, an estimate of MCE for a specific R is possible.}
}

\value{
estimates of MCE from bootstrap grouping replication plot along with the plot, if specified
}
\references{  }
\author{ Elizabeth Koehler and Sebastien Haneuse}
\note{ 
 This should be used cautiously when the Central Limit Theorem does not apply, or when the simulation-based estimate 
is not generated as an average of the replication results.  Note that the top axis presents the R that 
corresponds to the x-axis value of 1/sqrt(R) directly below.
}
\seealso{ \code{\link{mceCLT}}, \code{\link{mceBoot}}, \code{\link{RProj}}, \code{\link{groupbootMCE}} in package MCE.}
\examples{
percentbias<-rnorm(10000,15, 4)
groupbootMCE(percentbias, Rseq=c(100,500,2000), B=20, bb=10, targetR=c(15000, 10000), method="bootstrap", type="PB", truth=15, plot=TRUE)

%function(data, Rseq=c(100,500), B=20, bb=1,targetR=NA, method="bootstrap",truth=NULL,plot=FALSE, digits=2, xAxis=NULL, yAxis=NULL, zAxis=NULL, type="",...)
%{
%   R <- switch(type,
%              "mean"=length(data),
%              "PB"=length(data),
%              "SE"=length(data),
%              "RE"=nrow(data))

%  ##
%  if(is.null(targetR)) targetR <- R
%  ave1 <- matrix(NA,nrow=length(Rseq), ncol=bb)

%  for(k in 1:bb){
%    dataperm <- sample(1:R, R, replace=FALSE)
%    groups <- trunc(R/Rseq)
  
%    for(j in 1:length(Rseq)){
%      MCE <- rep(NA,groups[j])
%      data1 <- rep(NA,groups[j])
          
%  if(method=="clt")
%      for(i in 1:groups[j]){
%        data1[1:Rseq[j]] <- dataperm[(((i-1)*Rseq[j])+1):(i*Rseq[j])]
%        MCE[i] <- mceCLT(data[data1], type, truth)
%      }
        
%      if(method=="bootstrap")
%      for(i in 1:groups[j]){
%        data1[1:Rseq[j]] <- dataperm[(((i-1)*Rseq[j])+1):(i*Rseq[j])]
%        MCE[i] <- mceBoot(data[data1], B, type,truth)
%      }
%    ave1[j,k]<-mean(MCE)
%    }
%  }
%  Ymce <- ave1
%  Xmce <- rep(1 / sqrt(Rseq), each=bb)
%  beta1 <- lm(as.vector(t(Ymce)) ~ as.vector(Xmce) - 1)$coef
%  value <- as.numeric(beta1 / sqrt(targetR))   
% 
%  if(plot == TRUE)
%  {
%    ##
%    xRange <- c(0, ifelse(is.null(xAxis) == TRUE, 1/sqrt(min(Rseq)), max(xAxis))) 
%    yRange <- c(0, ifelse(is.null(yAxis) == TRUE, max(Ymce), max(yAxis)))
%    plot(xRange, yRange, axes=FALSE, type="n", xlab=expression(1 / sqrt(R)), ylab="Monte Carlo error")
%    ##
%    points(Xmce, as.vector(t(Ymce)), pch="X", cex=1.5)
%    if(is.null(xAxis) == TRUE) xAxis <- round(seq(from=0, to=xRange[2], length=5), digits=digits[1])
%    axis(1, at=xAxis)
%    if(is.null(yAxis) == TRUE) yAxis <- round(seq(from=0, to=yRange[2], length=5), digits=digits[1])
%    axis(2, at=yAxis)
%    if(is.null(zAxis) == TRUE) zAxis <- c(targetR, Rseq)
%    axis(3, at=1/sqrt(zAxis), labels=zAxis)
%    abline(0, beta1, col="red", lty=2)    
%    ##
%    segments(0, value,1/sqrt(targetR), value, col="blue", lwd=2)
%    segments(1/sqrt(targetR), 0, 1/sqrt(targetR), value, col="blue", lwd=2)
%    text(1/sqrt(targetR[1]) * 1.1, mean(yAxis[(length(yAxis)-c(0,1))])-seq(0,(length(value)-1)*.25*yAxis[2], by=.25*yAxis[2]), paste("MCE estimate at R=",sort(targetR, decreasing=FALSE),"=", round(sort(value, decreasing=TRUE), digits=digits[length(digits)]), sep=" "), pos=4)
%    box(bty="l")
%  }
%  return(value)
%  }
}

\keyword{univar  }

