`mceBoot` <-
function(data, B, type="", truth=NULL)
{
  ##
  R <- switch(type,
              "mean"=length(data),
              "PB"=length(data),
              "SE"=length(data),
              "RE"=nrow(data))

  ##
  subsets <- matrix(sample(1:R, R*B, replace=TRUE), ncol=B)
  bootRes <- rep(NA, B)
  for(b in 1:B)
  {
    if(type == "mean") bootRes[b] <- calcEX(data, subsets[,b])
    if(type == "PB")   bootRes[b] <- calcPB(data, subsets[,b], truth=truth)
    if(type == "SE")   bootRes[b] <- calcSE(data, subsets[,b])
    if(type == "RE")   bootRes[b] <- calcRE(data, subsets[,b])
  }
  value <- sd(bootRes)
  return(value)
}

