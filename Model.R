#' The functions in this file prepare the models
#' @param miss is the maximum percentage of missing values allowed in the feature. All features with more missing values will be dropped from the analysis process.
#' @param correlation is the maximum positive and negative correlation allowed in the predictor correlation matrix. All features with higher correlation values are dropped from the analysis process.
#' @param dataset is the code for the dataset for which this study will be performed. 1: is for CHSI dataset and 3: is for SWAN dataset.
#' @param seed is for result repetition.
#' @param sample_clustersize is the number of pieces the samples has to be fragmented for DMU.
#' @param sample_to_featureratio is the minimum sample to feature ratio allowed in each of the cluster. Any cluster with less than the minimum sample to feature ratio is dropped from the analysis process.
#' @param datatype determines whether the training data will have any sample with complete data or not. "Both" means that training data will have samples containing complete and incomplete data. "Incomplete" means that training data will have samples containing only incomplete data.
#'
#'
#' @export
fit=function(x = list(miss=10, correlation=0.52, dataset=1, seed=1, sample_clustersize = 8, sample_to_featureratio=2, datatype="incomplete"), technique= "DMU"){
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

  # Run the Models
  if(length(dataset$train)>1){
    # Create data split for DMU
    traindata = dataset$train
    testdata = dataset$test
    print(nrow(traindata))
    print(nrow(testdata))
    # Define Hyperpara to label the df
    hyperpara=data.frame(miss = missing_percent, corr=corr_range, cutter=ncol(traindata)-1, trainset=datatype, seed=seed, data_code=data_code, stringsAsFactors = F)
    # Run the Model
    if(technique == "reg"){
      data_prepare=list(list(NA, comp_set = testdata,incomplete = traindata))
      res=mreg_results(hyperpara, data_prep=data_prepare)
    }
    else if(technique == "mice"){
      data_prepare=list(list(NA, comp_set = testdata,incomplete = traindata))
      res=mmice_results(hyperpara, data_prep=data_prepare)
    }
    else{
      dfsplit=mdatasplit(traindata, cutratio=cutratio, cluster_number=sample_clustersize)
      data_prepare=list(list(dfsplit,comp_set = testdata, incomplete = traindata))
      res=mDMU(hyperpara, data_prep=data_prepare)
    }
  }else{print("Error: Check the training data")}
  fin_res=res[[1]]
  return(fin_res)
}
