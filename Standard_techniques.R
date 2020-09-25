#' The funcstions in this file focus determining the performance of the models using linear regression and predictive mean matching based MICE.


# Run the regression
reg_results=function(hyperpara, data_prep){
  reg_res=lapply(1:(nrow(hyperpara)-0), function(x) {
    miss=hyperpara[x,1]
    corr=hyperpara[x,2]
    cutter=hyperpara[x,3]
    train=hyperpara[x,4]
    seed=hyperpara[x,5]
    data_code=hyperpara[x,6]
    data=data_prep[[x]]
    if(length(data) == 3){
      if(train=="both" & (nrow(data[[2]])/ncol(data[[2]]))>9){
        set.seed(seed);index=sample(1:nrow(data[[2]]), floor(nrow(data[[2]])*5/10))
        testset=data[[2]][index,]
        train_comp=data[[2]][-index,]
        len_test_y=length(unique(testset$y))
        len_train_y=length(unique(train_comp$y))
        if(len_test_y>1 & len_train_y>1){
          inputdata=testset
          trainset=rbind(train_comp,data[[3]])
          reg=lm(y~., data = trainset)
          lm_reg=data.frame(beta_Reg=summary(reg)$coefficients[,1],
                            sd_lm=summary(reg)$coefficients[,2],
                            Variable=attributes(summary(reg)$coefficients)$dimnames[[1]],
                            stringsAsFactors = F)
          y_estimate=y_est_nointerac(betadf = lm_reg,xtestdf = inputdata, meanbeta = "beta_Reg")
          res=predict_metric(actualy = inputdata[,"y"], predictedy=y_estimate)
          res$Approach = "beta_Reg"; res$traindata = train; res$missing = miss;
          res$correlation = corr; res$sampleratio = cutter; res$nearzerovar = 10;
          res$seed = seed;res$split=0; res$file=data_code
          res
        }else{res=NULL}
      }else{res=NULL}
    }else{res=NULL}
  })
  return(reg_res)
}

# Run the Mice regression
mice_results=function(hyperpara, data_prep){
  mice_res=lapply(1:(nrow(hyperpara)), function(x) {
    miss = hyperpara[x,1]
    corr = hyperpara[x,2]
    cutter = hyperpara[x,3]
    train = hyperpara[x,4]
    seed = hyperpara[x,5]
    data_code = hyperpara[x,6]
    data = data_prep[[x]]
    if(length(data) == 3){
      if((nrow(data[[2]])/ncol(data[[2]]))>9){
        if(train=="both"){
          set.seed(seed);index=sample(1:nrow(data[[2]]), floor(nrow(data[[2]])*5/10))
          testset=data[[2]][index,]
          train_comp=data[[2]][-index,]
          len_test_y=length(unique(testset$y))
          len_train_y=length(unique(train_comp$y))
          if(len_test_y>1 & len_train_y>1){
            inputdata=testset
            trainset=rbind(train_comp,data[[3]])
            res=1
          }else{res=NULL}
        }
        else{
          ## Train has only incomplete data
          trainset=data[[3]]
          testset=data[[2]]
          res=1
        }
      }else{res=NULL}
    }else{res=NULL}
    if(!is.null(res)){
      inputdata=testset
      drop_yindex=which(colnames(trainset)=="y")
      imp=mice::mice(trainset[,c(-drop_yindex)], seed = 1, printFlag = F)
      mice_lm_res=lapply(1:5, function(imputer){
        comp_imp=mice::complete(imp,imputer)
        Var_model=names(comp_imp)
        comp_imp$y=trainset$y
        mice_res=lm(y~., data=comp_imp)
        res=summary(mice_res)$coefficients[,1]
        res
      })

      mice_lm=t(apply(data.frame(mice_lm_res),1,function(x) {c(beta_Imputed_Reg=mean(x, na.rm=TRUE),sd_mice=sd(x, na.rm = TRUE))}))
      mice_lm=data.frame(mice_lm, Variable=rownames(mice_lm), stringsAsFactors = F)
      mice_lm[is.na(mice_lm)]=0
      y_estimate=y_est_nointerac(betadf = mice_lm,xtestdf = inputdata, meanbeta = "beta_Imputed_Reg")
      res=predict_metric(actualy = inputdata[,"y"], predictedy=y_estimate)
      res$Approach="beta_Imputed_Reg"
      res$traindata=train; res$missing=miss; res$correlation=corr; res$sampleratio=cutter;
      res$nearzerovar=10; res$seed=seed;res$split=0; res$file=data_code
      res
    }
  })
  return(mice_res)
}

mreg_results=memoise::memoise(reg_results)
mmice_results=memoise::memoise(mice_results)
