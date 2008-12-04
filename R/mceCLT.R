`mceCLT` <-
function(data, type="", truth=NULL)
{
  ##
  R <- length(data)

  ##
  value <- NA
  if(type == "mean") value <- sd(data) / sqrt(R)
  if(type == "PB") value <- sd((data - truth) / truth * 100) / sqrt(R)
  return(value)
}

