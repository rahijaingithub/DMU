#' The funcstions in this file focus determining the performance of the models
#'
# Estimate the output
y_est_nointerac=function(betadf, xtestdf, meanbeta){
  #print(betadf)
  # Get Intercept
  Intercept=ifelse(any(betadf$Variable=="(Intercept)")==TRUE,
                   betadf[which(betadf$Variable=="(Intercept)"), meanbeta], 0)
  betadf=data.frame(betadf[, c("Variable", meanbeta)])

  #Multiply Beta with data
  multi=data.frame(lapply(betadf$Variable[which(betadf$Variable!="(Intercept)")], function(x) {
    a=betadf[which(betadf$Variable==x), meanbeta]
    b= xtestdf[,x]
    a*b
  }))
  multi=data.frame(do.call(cbind,multi))

  # Add Intercept
  multi=data.frame(cbind(multi,Intercept))

  #Get Y
  y_est=apply(multi, 1, sum)
  return(y_est)
}
rsquared=function(actualy, predictedy){
  mean_actual=mean(actualy, na.remove=T)
  rss=sum((actualy-predictedy)^2)
  ess=sum((mean_actual-predictedy)^2)
  tss=ess+rss
  rsquare=1-(rss/tss)
  output=c(rsquare,rss,ess)
  return(rsquare)
}
predict_metric=function(actualy, predictedy){
  Corr= cor(actualy,predictedy, method = "pearson")
  MSE= MLmetrics::MSE(actualy,predictedy)
  rsquare=rsquared(actualy = actualy, predictedy = predictedy)
  res=data.frame(Corr=Corr,MSE=MSE, rsquare=rsquare, stringsAsFactors = F)
  return(res)
}
