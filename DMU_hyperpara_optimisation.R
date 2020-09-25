#' The function in this file perform hyperparameter optimisation for DMU.

#' @export
hyper_para_optim = function(x = list(miss=10, correlation=0.52, dataset=1, seed=1, sample_clustersize = NA, sample_to_featureratio=2, datatype="incomplete"), traindata, testdata){
  #Decode
  missing_percent=x[[1]]
  corr_range=x[[2]]
  data_code=x[[3]]
  seed=x[[4]]
  sample_clustersize=x[[5]]
  cutratio=x[[6]]
  datatype=x[[7]]

  # Generate the dataset
  hyperpara = data.frame(miss= missing_percent, corr= corr_range, seed =1, data_code= data_code, stringsAsFactors = F)
  dataset = data_extract(hyperparameter = hyperpara)
  traindata = dataset$train
  testdata = dataset$test

  # Define Hyperpara to label the df
  hyperpara = data.frame(miss = missing_percent, corr=corr_range, cutter=ncol(traindata)-1, trainset=datatype, seed=seed, data_code=data_code, stringsAsFactors = F)

  maximum_clusters=nrow(traindata)/cutratio
  GA_FUNCT=function(bin){
    dec_code = GA::binary2decimal(bin)
    if(dec_code<(maximum_clusters+1) & dec_code>0){
      dfsplit=mdatasplit(traindata, cutratio=cutratio, cluster_number=dec_code)
      data_prepare=list(list(dfsplit,testdata,traindata))
      bay_res=mDMU(hyperpara, data_prep=data_prepare)
      if(length(bay_res[[1]]) != 0){
        if(!is.na(bay_res[[1]]$MSE)){res=-bay_res[[1]]$MSE}
        else{res=-10}}
      else{res=-10}
    }else(res=-10)
    return(Score=res)
  }
  nBits = length(GA::decimal2binary(maximum_clusters))
  set.seed(3)
  GA <<- GA::ga(type = "binary", fitness = GA_FUNCT, nBits = nBits, popSize = 10, maxiter = floor(maximum_clusters/20), maxFitness = 0.1*100000, run=10, pcrossover = 0.8, pmutation = 0.8, parallel = F, monitor = F)
  best_cluster = GA::binary2decimal(GA@solution[1,])
  return(best_cluster)
}
