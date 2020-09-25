#' The functions in this file prepare the dataset for the DMU modeling

#' # Remove the variables with missing data
#'
library(magrittr)
missing_data=function(inputdf, missingpercent){
  # Determine the missing data distribution across the dataset
  Missing_data=sapply(inputdf, function(x) 100*sum(is.na(x))/nrow(inputdf))
  #hist(Missing_data, xlab="Percentage Missing Data", ylab="Number of Input Variables", main="Histogram of Missing Data")
  # Remove variables with more than "missingpercent"
  var=attributes(Missing_data)$names[Missing_data>missingpercent]
  inputdf[,var]=NULL
  return(inputdf)
}
# Remove the variables with Low variation
lowvar=function(inputdf, minfreq=100/10){
  al = unlist(apply(inputdf,2,function(x) caret::nearZeroVar(x, freqCut = minfreq, saveMetrics = T)[4]))
  drop_var_2 = names(al)[which(al==TRUE)]
  drop_var_2 = stringr::str_replace(drop_var_2,".nzv","")
  #print("newdf?")
  newdf=inputdf[,setdiff(names(inputdf),drop_var_2)]
  #print("newdf found")
  return(newdf)
}
# Remove the variables with High correlation
Correl=function(inputdf, cutoff){
  dfcor=cor(inputdf, use = "pairwise.complete.obs")
  highcor=apply(dfcor,1,function(x) length(which(x>cutoff | x< (-1*cutoff))))
  var=attributes(highcor)$names
  newdf=inputdf[,var]
  while (max(highcor>1)) {
    sort=order(-highcor)
    # drop the variable
    var=attributes(highcor)$names[sort[-1]]
    #print(var)
    newdf=inputdf[,var]
    dfcor=cor(newdf, use = "pairwise.complete.obs")
    highcor=apply(dfcor,1,function(x) length(which(x>cutoff | x< (-1*cutoff))))
    #print(c("list",highcor))
  }
  #print("exit")
  return(newdf)
}
# Remove the Outliers
out_rem=function(inputdf){
  Meanprofile=apply(inputdf, 2, function(x) mean(x, na.rm=T))
  sdprofile=apply(inputdf, 2, function(x) sd(x, na.rm=T))
  m_plus_s=Meanprofile+(4*sdprofile)
  m_minus_s=Meanprofile-(4*sdprofile)
  for(x in names(inputdf)){
    len=length(inputdf[which(inputdf[,x]>m_plus_s[x] | inputdf[,x]<m_minus_s[x] ),x])
    #print(len)
    if(len>0){
      inputdf=inputdf[-which(inputdf[,x]>m_plus_s[x] | inputdf[,x]<m_minus_s[x]),]
    }
    #print(nrow(inputdf))
  }
  return(inputdf)
}

# Data Preparation Function
#' @export
df_process=function(inputdf,outcomedf,missingper=70,minvar=100/10,corr_cut=0.1, outcomevar=c("Unhealthy_Days")){
  # Remove the independent variables with missing data and near zero variable frequency
  input_miss=missing_data(inputdf=inputdf, missingpercent=missingper)
  # str(input_miss)
  #print("Found miss_data")
  input_nozero=lowvar(inputdf=input_miss, minfreq=minvar)
  # str(input_nozero)

  # Remove variables with high correlations
  input_nocorrel=Correl(inputdf=input_nozero, cutoff=corr_cut)
  # str(input_nocorrel)
  dataset=data.frame(apply(input_nocorrel,2, scales::rescale))
  # str(dataset)
  dataset$y=outcomedf[,outcomevar]
  final_df=dataset[which(!is.na(dataset$y)),]

  comp_set=final_df[complete.cases(final_df),]
  incomplete=final_df[!complete.cases(final_df),]

  # Remove the empty columns, low variance columns and dataset with no "y" from final dataset
  incomplete_pre=incomplete %>% janitor::remove_empty_cols()
  incomplete=lowvar(inputdf=incomplete_pre, minfreq=minvar)
  comp_set_pre=comp_set[,names(incomplete)]
  comp_set_pre=lowvar(inputdf=comp_set_pre, minfreq=minvar)

  if(ncol(comp_set_pre)>1){
    if(names(comp_set_pre)[ncol(comp_set_pre)]=="y"){
      comp_set=comp_set_pre
      incomplete=incomplete[,names(comp_set)]
      return(list(train=incomplete, test=comp_set, cormat= NA))
    }else{return(list(train=NULL, test=NULL, cormat= NA))}
  }
  else{return(list(train=NULL, test=NULL, cormat= NA))}

}
#' @export
data_extract=function(hyperparameter=data.frame(missing=30, correlation=0.52, seed =1, dataset=1), seed=1){
  #print(seed)
  data_prep=lapply(1:(nrow(hyperparameter)), function(x) {
    miss=hyperparameter[x,1]
    corr=hyperparameter[x,2]
    if(is.na(seed)){seed=hyperparameter[x,3]}

    # Get the dataset
    {
      data_code=hyperparameter[x,4]
      in_data=  paste0("https://raw.githubusercontent.com/rahijaingithub/Dynamic-Model-Updating-DMU-/master/Dataset/in_d",data_code,".csv")
      out_data= paste0("https://raw.githubusercontent.com/rahijaingithub/Dynamic-Model-Updating-DMU-/master/Dataset/out_d",data_code,".csv")
      Outdata=read.csv(out_data)
      Inputdata=read.csv(in_data)
      Outdata=data.frame(out=Outdata[,2])
      Inputdata=Inputdata[,-1]
    }
    dp=mdf_process(missingper=miss, corr_cut=corr, inputdf=Inputdata,outcomedf=Outdata,minvar=100/10, outcomevar=c("out"))
    dp
  })
  # str(data_prep)
  return(data_prep[[1]])
}

# Memoisation of functions
mdf_process = memoise::memoise(df_process)
mdata_extract = memoise::memoise(data_extract)
