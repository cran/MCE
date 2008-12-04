`mceProj` <-
function(data, Rseq, B, bb=1,targetR=NULL, type="", truth=NULL, plot=FALSE, digits=2, xAxis=NULL, yAxis=NULL, zAxis=NULL, ...)
{
  ##
  R <- switch(type,
              "mean"=length(data),
              "PB"=length(data),
              "SE"=length(data),
              "RE"=nrow(data))

  ##
  if(is.null(targetR)) targetR <- R
  
  ##
  Ymce  <- matrix(NA, nrow=bb, ncol=length(Rseq))
  Xmce  <- matrix(NA, nrow=bb, ncol=length(Rseq))
  for(k in 1:bb){
  bootRes <- matrix(NA, nrow=B, ncol=length(Rseq))
  for(i in 1:length(Rseq))
  {
    subsets <- matrix(sample(1:R, Rseq[i]*B, replace=TRUE), ncol=B)
    for(b in 1:B)
    {
      if(type == "mean") bootRes[b,i] <- calcEX(data, subsets[,b])
      if(type == "PB")   bootRes[b,i] <- calcPB(data, subsets[,b], truth=truth)
      if(type == "SE")   bootRes[b,i] <- calcSE(data, subsets[,b])
      if(type == "RE")   bootRes[b,i] <- calcRE(data, subsets[,b])
    }
  }
  ##
  Ymce[k,]  <- apply(bootRes, 2, sd)
  Xmce[k,]  <- 1 / sqrt(Rseq)
  ##
}
  beta1 <- lm(as.vector(Ymce) ~ as.vector(Xmce) - 1)$coef
  value <- as.numeric(beta1 / sqrt(targetR))
  if(plot == TRUE)
  {
    ##
    xRange <- c(0, ifelse(is.null(xAxis) == TRUE, 1/sqrt(min(Rseq)), max(xAxis)))
    yRange <- c(0, ifelse(is.null(yAxis) == TRUE, max(Ymce), max(yAxis)))
    plot(xRange, yRange, axes=FALSE, type="n", xlab=expression(1 / sqrt(R)), ylab="Monte Carlo error", ...)
    ##
    points(Xmce, Ymce, pch="X", cex=1.5)
    if(is.null(xAxis) == TRUE) xAxis <- round(seq(from=0, to=xRange[2], length=5), digits=digits[1])
    axis(1, at=xAxis)
    if(is.null(yAxis) == TRUE) yAxis <- round(seq(from=0, to=yRange[2], length=5), digits=digits[1])
    axis(2, at=yAxis)
    if(is.null(zAxis) == TRUE) zAxis <- c(targetR, Rseq)
    axis(3, at=1/sqrt(zAxis), labels=zAxis)
    abline(0, beta1, col="red", lty=2)    
    ##
    segments(0, value, 1/sqrt(targetR), value, col="blue", lwd=2)
    segments(1/sqrt(targetR), value, 1/sqrt(targetR), 0, col="blue", lwd=2)
    text(1/sqrt(targetR[1]) * 1.1, mean(yAxis[(length(yAxis)-c(0,1))])-seq(0,(length(value)-1)*.25*yAxis[2], by=.25*yAxis[2]), paste("MCE estimate =", round(sort(value, decreasing=TRUE), digits=digits[length(digits)]), sep=" "), pos=4)
    box(bty="l")
  }
  return(value)
}

