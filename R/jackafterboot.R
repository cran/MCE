`jackafterboot` <-
function(data, B=100,type="", truth=NULL){
  R <- switch(type,
              "mean"=length(data),
              "PB"=length(data),
              "SE"=length(data),
              "RE"=nrow(data))

     x<-seq(1:B)
     est<-rep(NA, length(data))
     
     subsets <- matrix(sample(1:length(data), length(data)*B, replace=TRUE), nrow=B, byrow=TRUE)
      
     for(i in 1:length(data)){

          pts<-which(subsets==i, arr.ind=TRUE)
          exclude<-unique(pts[,1])
          include<-setdiff(x,exclude)
          jack<-subsets[include,]
          ave<-rep(NA, length(include))
          bootRes <- rep(NA, length(include))
          for(j in 1:length(include)){
                if(type == "mean") bootRes[j] <- calcEX(data, subsets[include[j],])
                if(type == "PB")   bootRes[j] <- calcPB(data, subsets[include[j],], truth=truth)
                if(type == "SE")   bootRes[j] <- calcSE(data, subsets[include[j],])
                if(type == "RE")   bootRes[j] <- calcRE(data, subsets[include[j],])
          }
          est[i]<-sd(bootRes)
     }
     MCEjack<-sqrt((length(data)-1)*var(est))
     return(MCEjack)
}

