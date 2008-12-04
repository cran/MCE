`groupbootR` <-
function(data, Rseq=c(100,500), B=20, bb=1,targetMCE=NA, method="bootstrap",truth=NULL,plot=FALSE, digits=2, xAxis=NULL, yAxis=NULL, zAxis=NULL, type="",...)
{
   R <- switch(type,
              "mean"=length(data),
              "PB"=length(data),
              "SE"=length(data),
              "RE"=nrow(data))

  ##

  ave1 <- matrix(NA,nrow=length(Rseq), ncol=bb)

  for(k in 1:bb){
    dataperm <- sample(1:R, R, replace=FALSE)
    groups <- trunc(R/Rseq)
  
    for(j in 1:length(Rseq)){
      MCE <- rep(NA,groups[j])
      data1 <- rep(NA,groups[j])
          
  if(method=="clt")
      for(i in 1:groups[j]){
        data1[1:Rseq[j]] <- dataperm[(((i-1)*Rseq[j])+1):(i*Rseq[j])]
        MCE[i] <- mceCLT(data[data1], type, truth)
      }
        
      if(method=="bootstrap")
      for(i in 1:groups[j]){
        data1[1:Rseq[j]] <- dataperm[(((i-1)*Rseq[j])+1):(i*Rseq[j])]
        MCE[i] <- mceBoot(data[data1], B, type,truth)
      }
    ave1[j,k]<-mean(MCE)
    }
  }
  Ymce <- ave1
  Xmce <- rep(1 / sqrt(Rseq), each=bb)                
  beta1 <- lm(as.vector(t(Ymce)) ~ as.vector(Xmce) - 1)$coef
  value <- round(as.numeric((targetMCE/beta1)^(-2))) 

  if(plot == TRUE)
  {
    ##
    xRange <- c(0, ifelse(is.null(xAxis) == TRUE, 1/sqrt(min(Rseq)), max(xAxis))) 
    yRange <- c(0, ifelse(is.null(yAxis) == TRUE, max(Ymce), max(yAxis)))
    plot(xRange, yRange, axes=FALSE, type="n", xlab=expression(1 / sqrt(R)), ylab="Monte Carlo error")
    ##
    points(Xmce, as.vector(t(Ymce)), pch="X", cex=1.5)
    if(is.null(xAxis) == TRUE) xAxis <- round(seq(from=0, to=xRange[2], length=5), digits=digits[1])
    axis(1, at=xAxis)
    if(is.null(yAxis) == TRUE) yAxis <- round(seq(from=0, to=yRange[2], length=5), digits=digits[1])
    axis(2, at=yAxis)
    if(is.null(zAxis) == TRUE) zAxis <- c(value, Rseq)
    axis(3, at=1/sqrt(zAxis), labels=zAxis)
    abline(0, beta1, col="red", lty=2)    
    ##
    segments(0, targetMCE,1/sqrt(value), targetMCE, col="blue", lwd=2)
    segments(1/sqrt(value), 0, 1/sqrt(value), targetMCE, col="blue", lwd=2)
    text(1/sqrt(value[1]) * 1.1, mean(yAxis[(length(yAxis)-c(0,1))])-seq(0,(length(value)-1)*.25*yAxis[2], by=.25*yAxis[2]), paste("R estimate at MCE=",sort(targetMCE, decreasing=TRUE),"=", round(sort(value, decreasing=FALSE), digits=digits[length(digits)]), sep=" "), pos=4)
    box(bty="l")
  }
  return(value)
}

