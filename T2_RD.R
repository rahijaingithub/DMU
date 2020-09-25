# For Real Dataset

# Load packages
library(pacman)
pacman::p_load(readxl,foreach, doParallel, MCMCpack, mice, caret, stringr, MLmetrics, reshape2, ggplot2, janitor, Hmisc, 
               scales, psych, purrr, fGarch, mosaic, copula, gtools, GA, plyr, memoise, ParBayesianOptimization)

# Define the Functions:
{
  ## Prepare the data
  # Remove the variables with missing data, low variation and high correlation
  {
    missing_data=function(inputdf, missingpercent){
      # Determine the missing data distribution across the dataset
      Missing_data=sapply(inputdf, function(x) 100*sum(is.na(x))/nrow(inputdf))
      #hist(Missing_data, xlab="Percentage Missing Data", ylab="Number of Input Variables", main="Histogram of Missing Data")
      # Remove variables with more than "missingpercent"
      var=attributes(Missing_data)$names[Missing_data>missingpercent]
      inputdf[,var]=NULL
      return(inputdf)
    }
    lowvar=function(inputdf, minfreq=100/10){
      al= unlist(apply(inputdf,2,function(x) nearZeroVar(x, freqCut = minfreq, saveMetrics = T)[4]))
      drop_var_2=names(al)[which(al==TRUE)]
      drop_var_2=str_replace(drop_var_2,".nzv","")
      #print("newdf?")
      newdf=inputdf[,setdiff(names(inputdf),drop_var_2)]
      #print("newdf found")
      return(newdf)
    }
    ## To keep maximum variables
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
  # Cut the data for the proposed method
  datasplit=function(inputdf, cutratio=5, cluster_number=NA){
    #print(c("CutRatio", cutratio))
    norm_set=inputdf
    code_set=inputdf
    code_set[!(is.na(code_set))]<-1
    code_set[is.na(code_set)]<-0
    #print(cluster_number)
    ## Make rowcluster
    hc=hclust(dist(code_set),method="ward.D2")
    ### Optimize the number of cuts
    cut=floor(nrow(code_set)/cutratio)
    # scan for the cuts which give variable to samplesize ratio of atleast 5
    if(is.na(cluster_number)){
      cutlist=lapply(1:cut, function(no_of_clusters) {
        x = no_of_clusters
        #print(x)
        fit=cutree(hc, x)
        print(c("no_of_clusters: ", x))
        # Those cuts are processed which meet minimum cutratio criteria
        if(min(table(fit))<cutratio){out=NULL}
        else{
          print(min(table(fit)))
          norm_set$fit=fit
          code_set$fit=fit
          store=list()
          # For each cut size
          ini_list=lapply(1:x, function(i){
            
            clus_len=length(fit[which(fit==i)])
            a=code_set[which(code_set$fit==i),1:(ncol(code_set)-2)]
            
            a = a[colSums(a)==nrow(a)]
            constant_y=unique(norm_set[which(norm_set$fit==i),"y"])
            #print((nrow(a)/ncol(a)))
            # Make sure only those cuts are retained which does not have constant "y" and sufficient data
            if(ncol(a)>0 & (nrow(a)/ncol(a)) >= cutratio & length(constant_y)>1){
              tempdf=norm_set[which(norm_set$fit==i),c(names(a),"y")]
              tempdf=tempdf[lapply(tempdf, function(x) length(unique(x)))>1]
              if(ncol(tempdf)>1){store[[i]]=tempdf}
              else{store[[i]]=NULL}
            }else{store[[i]]=NULL}
            store=store[lapply(store,length)>0]
          })
          # print("Initial")
          # str(ini_list)
          ini_list=ini_list[lapply(ini_list,length)>0]
          out=ini_list
          out
        }
        out=out[lapply(out,length)>0]
        out
      })
    }
    else{
      cutlist=lapply(cluster_number, function(no_of_clusters) {
        x = no_of_clusters
        #print(x)
        if(cluster_number>cut | cluster_number<1){out=NULL}
        else{
          fit=cutree(hc, x)
          #print(fit)
          
          # Those cuts are processed which meet minimum cutratio criteria
          ## Step 1: Eliminate the cluster groups containing any cluster whose row to column ratio is less than cut ratio
          if(min(table(fit))<cutratio){out=NULL}
          else{
            #print(min(table(fit)))
            #print(table(fit))
            norm_set$fit=fit
            code_set$fit=fit
            #store=list()
            # For each cut size or cluster, find the variables with complete rows
            ini_list=lapply(1:x, function(i){
              #print(c("Cluster_name: ",i))
              numb_rows_in_clus = length(fit[which(fit==i)])
              #print(c(x, clus_len))
              
              a=code_set[which(code_set$fit==i),1:(ncol(code_set)-2)]
              #print(colSums(a))
              a = a[colSums(a)==nrow(a)]
              constant_y=unique(norm_set[which(norm_set$fit==i),"y"])
              #print((nrow(a)/ncol(a)))
              # Step 2: Make sure only those cuts are retained which does not have constant "y" and sufficient data
              if(ncol(a) > 0 & (nrow(a)/ncol(a)) >= cutratio & length(constant_y) > 1){
                tempdf=norm_set[which(norm_set$fit==i),c(names(a),"y")]
                #str(tempdf)
                tempdf=tempdf[lapply(tempdf, function(x) length(unique(x)))>1]
                #str(tempdf)
                if(ncol(tempdf)>1){store=tempdf}
                else{store=NULL}
              }else{store=NULL}
              #str(store)
              store
            })
            # print("Initial")
            #str(ini_list)
            ini_list=ini_list[lapply(ini_list,length)>0]
            #str(ini_list)
            out=ini_list
            out
          }
          out=out[lapply(out,length)>0]
          #str(out)
          out
        }
        
      })
      #cutlist=list(cutlist)
    }
    #str(cutlist)
    dfcut=cutlist[lapply(cutlist,length)>0]
    return(dfcut)
  }
  
  # Data Preparation Function
  df_process=function(inputdf,outcomedf,missingper=70,minvar=100/10,corr_cut=0.1, outcomevar=c("Unhealthy_Days")){
    # Remove the independent variables with missing data and near zero variable frequency
    input_miss=missing_data(inputdf=inputdf, missingpercent=missingper)
    #str(input_miss)
    #print("Found miss_data")
    input_nozero=lowvar(inputdf=input_miss, minfreq=minvar)
    #str(input_nozero)
    # Remove variables with high correlations
    input_nocorrel=Correl(inputdf=input_nozero, cutoff=corr_cut)
    #str(input_nocorrel)
    dataset=data.frame(apply(input_nocorrel,2, rescale))
    #str(dataset)
    dataset$y=outcomedf[,outcomevar]
    final_df=dataset[which(!is.na(dataset$y)),]
    
    comp_set=final_df[complete.cases(final_df),]
    incomplete=final_df[!complete.cases(final_df),]
    
    # Remove the empty columns, low variance columns and dataset with no "y" from final dataset
    incomplete_pre=incomplete %>% remove_empty_cols()
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
  
  data_extract=function(hyperparameter, seed=1){
    #print(seed)
    data_prep=lapply(1:(nrow(hyperparameter)), function(x) {
      #print(c("Para",x))
      miss=hyperparameter[x,1]
      corr=hyperparameter[x,2]
      if(is.na(seed)){seed=hyperparameter[x,3]}
      
      # Get the dataset
      {
        data_code=hyperparameter[x,4]
        #print(data_code)
        in_data=  paste0("~/in_d",data_code,".csv")
        out_data= paste0("~/out_d",data_code,".csv")
        #in_data=  paste0("D:/in_d",data_code,".csv")
        #out_data= paste0("D:/out_d",data_code,".csv")
        #in_data=  paste0("/gpfs/fs0/scratch/w/wxu/rahijain/in_d",data_code,".csv")
        #out_data= paste0("/gpfs/fs0/scratch/w/wxu/rahijain/out_d",data_code,".csv")
        Outdata=read.csv(out_data)
        Inputdata=read.csv(in_data)
        Outdata=data.frame(out=Outdata[,2])
        Inputdata=Inputdata[,-1]
        #print(ncol(Inputdata))
      }
      dp=mdf_process(missingper=miss, corr_cut=corr, inputdf=Inputdata,outcomedf=Outdata,minvar=100/10, outcomevar=c("out"))
      dp
    })
    return(data_prep[[1]])
  }
  
  # Estimate the output
  y_est_nointerac=function(betadf, xtestdf, meanbeta){
    #print(betadf)
    # Get Intercept
    Intercept=ifelse(any(betadf$Variable=="(Intercept)")==TRUE, 
                     betadf[which(betadf$Variable=="(Intercept)"), meanbeta], 0)
    #print(Intercept)
    
    betadf=data.frame(betadf[, c("Variable", meanbeta)])
    #Multiply Beta with data
    multi=data.frame(lapply(betadf$Variable[which(betadf$Variable!="(Intercept)")], function(x) {
      a=betadf[which(betadf$Variable==x), meanbeta]
      b= xtestdf[,x]
      #str(a*b)
      a*b
    }))
    multi=data.frame(do.call(cbind,multi))
    #str(multi)
    # Add Intercept
    multi=data.frame(cbind(multi,Intercept))
    #str(multi)
    #Get Y
    y_est=apply(multi, 1, sum)
    #print(y_est)
    return(y_est)
  }
  rsquared=function(actualy, predictedy){
    mean_actual=mean(actualy, na.remove=T)
    rss=sum((actualy-predictedy)^2)
    ess=sum((mean_actual-predictedy)^2)
    tss=ess+rss
    #print(c("rss", rss, "tss", tss))
    #plot(actualy,predictedy)
    rsquare=1-(rss/tss)
    output=c(rsquare,rss,ess)
    return(rsquare)
  }
  predict_metric=function(actualy=validationdf[,outname[x]], predictedy=y_estimate){
    Corr= cor(actualy,predictedy, method = "pearson")
    MSE= MSE(actualy,predictedy)
    rsquare=rsquared(actualy = actualy, predictedy = predictedy)
    res=data.frame(Corr=Corr,MSE=MSE, rsquare=rsquare, stringsAsFactors = F)
    return(res)
  }
  
  # Run the regression
  reg_results=function(hyperpara, data_prep){
    reg_res=lapply(1:(nrow(hyperpara)-0), function(x) {
      #reg_res=foreach::foreach(x=1:(nrow(hyperpara)-600),.packages = c("caret", "stringr", "janitor", "MLmetrics", "reshape2"))%dopar%{
      #print(c("reg", x))
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
  
  # Run the Bayesian
  bay_reg_mcmc=function(inputdf, outvar="y", df_beta, run=1){
    #print("bay enter")
    #str(inputdf)
    # Get df for the step and get predictor name
    Bay_Trial_Train=inputdf
    #print(c("outvar", outvar))
    #print(names(Bay_Trial_Train))
    Var_model=names(Bay_Trial_Train)[which(names(Bay_Trial_Train)!=outvar)]
    #print(Var_model)
    #Prepare the formula
    inputvar=paste0(Var_model,collapse = "+")
    #print(inputvar)
    equation=as.formula(paste0(c(outvar,inputvar), collapse = "~"))
    #print(equation)
    # Prepare the input values
    prior_beta=subset(df_beta, Variable %in% c("(Intercept)",Var_model))$Prior_mu
    prior_beta_precision=1/subset(df_beta, Variable %in% c("(Intercept)",Var_model))$Prior_sd^2
    c0=5; d0=50
    # Run the model and get the results
    #print(prior_beta)#;print(prior_beta_precision);print(equation); print(Bay_Trial_Train)
    #str(Bay_Trial_Train)
    model=MCMCregress(equation, data=Bay_Trial_Train, mcmc=5000, seed=1,b0=prior_beta, B0=prior_beta_precision, c0=c0, d0=d0)
    #print("success")
    out=summary(model)$statistics
    # Store the results and update the prior values
    out_df=data.frame(Variable=c("(Intercept)", Var_model),
                      Prior_mu=out[1:(nrow(out)-1),"Mean"], 
                      Prior_sd=out[1:(nrow(out)-1),"SD"],
                      posterior_mean=out[1:(nrow(out)-1),"Mean"], 
                      posterior_sd=out[1:(nrow(out)-1),"SD"], 
                      run=run, stringsAsFactors = F)
    #print(out_df)
    update_var=c("Prior_mu", "Prior_sd","posterior_mean","posterior_sd","run")
    df_beta[match(out_df$Variable,df_beta$Variable),update_var]=out_df[,update_var]
    return(df_beta)
  }
  bay_results=function(hyperpara, data_prep){
    #bay_res=lapply(1:(nrow(hyperpara)-640), function(x) {
    bay_res=foreach::foreach(x=1:(nrow(hyperpara)),.packages = c("caret", "stringr", "janitor", "MLmetrics", "reshape2"))%dopar%{
      #print(x)
      miss=hyperpara[x,1]
      corr=hyperpara[x,2]
      cutter=hyperpara[x,3]
      train=hyperpara[x,4]
      seed=hyperpara[x,5]
      data_code=hyperpara[x,6]
      data=data_prep[[x]]
      #str(data)
      if(length(data) == 3){
        if((nrow(data[[2]])/ncol(data[[2]]))>9 & length(data[[1]])>0){
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
        datacut=data[[1]]
        #str(datacut)
        #print(c("number of clusters: ", length(datacut)))
        
        bay_resul=lapply(1:length(datacut), function(round){
          #bay_res=foreach::foreach(round=1:length(datacut),.packages = c("MCMCpack", "mice", "caret", "stringr", "MLmetrics", "reshape2", "janitor"))%dopar%{
          #str(datacut[[1]])
          dflist=datacut[[round]]#[[1]]
          cut=length(dflist)
          #print(round)
          #print(dflist)
          #str(dflist)
          #Sys.sleep(2)
          if(train=="both"){
            dflist[[(cut+1)]]=train_comp
          }
          #str(dflist)
          var_model=names(trainset)[which(names(trainset)!="y")]
          Variable=union("(Intercept)", var_model)
          df_beta=data.frame(Variable=Variable, Prior_mu=0, Prior_sd=1e1, posterior_mean=NA, posterior_sd=NA, run=0, stringsAsFactors = F)
          #print(length(dflist))
          # Run Bayesian
          for(l in 1:length(dflist)){
            miss_data=data.frame(dflist[[l]])
            #str(miss_data)
            miss_data_var=intersect(names(trainset), names(miss_data))
            #print(c(l,miss_data_var))
            df_beta=bay_reg_mcmc(inputdf= miss_data[,miss_data_var], outvar = "y", df_beta = df_beta, run=l)
            #print(l)
          }
          bay_reg = data.frame(Variable=Variable,beta_Prop=df_beta$posterior_mean, sd_prop=df_beta$posterior_sd,  stringsAsFactors = F)
          bay_reg[is.na(bay_reg)]=0
          
          numb_varb=bay_reg$beta_Prop[1:5]
          names(numb_varb)=bay_reg$Variable[1:5]
          #print(numb_varb)
          
          y_estimate=y_est_nointerac(betadf = bay_reg,xtestdf = inputdata, meanbeta = "beta_Prop")
          #print(c(x,"Done so far"))
          bay_r=predict_metric(actualy = inputdata[,"y"], predictedy=y_estimate)
          bay_r$Approach="beta_Prop"
          bay_r$traindata=train; bay_r$missing=miss; bay_r$correlation=corr; 
          bay_r$sampleratio=cutter;bay_r$nearzerovar=10; bay_r$seed=seed;
          bay_r$split=length(dflist)
          #print(bay_r)
          bay_r
        })         
        res=do.call(rbind, bay_resul)
        res$file=data_code
        #print(res)
        res
      }
    }#)
    
    return(bay_res)
  }
  
  # Run the Mice regression
  mice_results=function(hyperpara, data_prep){
    mice_res=lapply(1:(nrow(hyperpara)), function(x) {
      #mice_res=foreach::foreach(x=1:(nrow(hyperpara)-600),.packages = c("mice","caret", "stringr", "janitor", "MLmetrics", "reshape2"))%dopar%{
      #print(hyperpara[x,])
      #print(x)
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
        #str(trainset)
        drop_yindex=which(colnames(trainset)=="y")
        #trainset<<-trainset[,c(-drop_yindex)]
        imp=mice(trainset[,c(-drop_yindex)], seed = 1, printFlag = F)
        #print(c(x, "success"))
        mice_lm_res=lapply(1:5, function(imputer){
          comp_imp=mice::complete(imp,imputer)
          # Create Interaction Variables and data
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
        #print(describe(res))
        res
      }
    })
    return(mice_res)
  }
}


# Memoise the functions
mdatasplit=memoise::memoise(datasplit)
mdf_process = memoise::memoise(df_process)
mdata_extract = memoise::memoise(data_extract)
mreg_results=memoise::memoise(reg_results)
mbay_results=memoise::memoise(bay_results)
mmice_results=memoise::memoise(mice_results)

fit=function(x){
  #print(x)
  #Decode
  missing_percent=x[[1]]
  corr_range=x[[2]]
  data_code=x[[3]]
  seed=x[[4]]
  run_GA=x[[5]]
  cutratio=x[[6]]
  datatype=x[[7]]
  
  # Generate the dataset
  hyperpara = data.frame(miss= missing_percent, corr= corr_range, seed =1, data_code= data_code, stringsAsFactors = F)
  dataset = data_extract(hyperparameter = hyperpara)
  
  # if(datatype == "incomplete"){testsize= 1000}else{testsize= 1050}
  # dataset=mdata_sim(varnum = x1, setting="Correlation", samplesize=(x1*x4)+testsize, testsize= testsize, effect = "Mar", Missingdata = T, 
  #                   maxmiss_per=x2, metric = T, corr_seed=x3, miss_seed=x3)
  
  # Run the Models
  if(length(dataset$train)>1){
    # Create data split for the Bayesian
    incomplete=dataset$train
    comp_set=dataset$test
    print(nrow(incomplete))
    print(nrow(comp_set))
    # Define Hyperpara to label the df
    hyperpara=data.frame(miss = missing_percent, corr=corr_range, cutter=ncol(incomplete)-1, trainset=datatype, seed=seed, data_code=data_code, stringsAsFactors = F)
    
    # Run the regression and store the results
    {
      data_prepare=list(list(NA,comp_set,incomplete))
      reg_res=mreg_results(hyperpara, data_prep=data_prepare)
      #print(reg_res)
      #return(reg_res)
    }
    # Run MICE regression and store the results
    {
      data_prepare=list(list(NA,comp_set,incomplete))
      mice_res=mmice_results(hyperpara, data_prep=data_prepare)
      #print(mice_res)
    }    
    # Run Bayreg and store the results
    # This is for fast modelng using genetic algorithm
    if(run_GA){
      #print("Enter")
      maximum_clusters=nrow(incomplete)/cutratio
      GA_FUNCT=function(bin){
        #print(bin)
        dec_code=binary2decimal(bin)
        #print(c("Dec_code:", dec_code, " max_cluster:", maximum_clusters))
        if(dec_code<(maximum_clusters+1) & dec_code>0){
          #res=fit(c(paralist[dec_code,]))
          #print("Enter")
          dfsplit=mdatasplit(incomplete, cutratio=cutratio, cluster_number=dec_code)
          data_prepare=list(list(dfsplit,comp_set,incomplete))
          bay_res=mbay_results(hyperpara, data_prep=data_prepare)
          #print(bay_res)
          #print(length(bay_res[[1]]))
          if(length(bay_res[[1]]) != 0){ 
            if(!is.na(bay_res[[1]]$MSE)){res=-bay_res[[1]]$MSE} #; print(res)
            else{res=-10}}
          else{res=-10}
        }else(res=-10)
        return(Score=res)
      }
      nBits=length(decimal2binary(maximum_clusters))
      set.seed(3)
      #print(floor(maximum_clusters/15))
      GA <<- ga(type = "binary", fitness = GA_FUNCT, nBits = nBits, popSize = 10, maxiter = floor(maximum_clusters/20), maxFitness = 0.1*100000, run=10, pcrossover = 0.8, 
                pmutation = 0.8, parallel = F, monitor = F)
      best_cluster=binary2decimal(GA@solution[1,])
      #print(GA@solution)
      #print(decimal2binary(75))
      #print(best_cluster)
      dfsplit=mdatasplit(incomplete, cutratio=cutratio, cluster_number=best_cluster)
      data_prepare=list(list(dfsplit,comp_set,incomplete))
      bay_res=mbay_results(hyperpara, data_prep=data_prepare)
      #print(bay_res)
    } 
    else{
      {
        dfsplit = mdatasplit(incomplete, cutratio=cutratio, cluster_number = NA)
        data_prepare = list(list(dfsplit,comp_set,incomplete))
        bay_res = mbay_results(hyperpara, data_prep=data_prepare)
        #print(bay_res)
      }
    }
    
    if(datatype == "both"){
      if(!is.null(reg_res[[1]]) & !is.null(bay_res[[1]]) & !is.null(mice_res[[1]])){
        reg_diff=reg_res[[1]]$MSE - min(bay_res[[1]]$MSE)
        mice_diff=mice_res[[1]]$MSE - min(bay_res[[1]]$MSE)
        f=ifelse(min(reg_diff,mice_diff)>0, reg_diff*mice_diff, -0.0000009)
      }
      else{f=-0.0000004}
    }
    else{
      #For incomplete
      if(!is.null(bay_res[[1]]) & !is.null(mice_res[[1]])){
        mice_diff=mice_res[[1]]$MSE-min(bay_res[[1]]$MSE)
        #print(mice_diff)
        f=ifelse(min(mice_diff)>0, mice_diff, -0.0000009)
      }
      else{f=-0.0000004}
    }
  }else{f=-0.0000001}
  
  #print(f)
  #getres <<- f
  if(datatype == "both"){
    flist=list(reg_res[[1]], bay_res[[1]][which(bay_res[[1]]$MSE==min(bay_res[[1]]$MSE)),], mice_res[[1]])
  }
  else{flist=list(bay_res[[1]][which(bay_res[[1]]$MSE==min(bay_res[[1]]$MSE)),], mice_res[[1]])}
  #print(mice_res)
  
  flist=do.call(rbind, flist)
  #print(flist)
  return(flist)
}
mfit=memoise::memoise(fit)
summarizer_bay=function(x, index=1){
  #print(x)
  res_sum=x %>% dplyr::group_by(Approach, traindata) %>%
    dplyr::summarise(missing_per=mean(as.numeric(missing), na.rm = T),
                     correlation=mean(as.numeric(correlation), na.rm = T),
                     predictor=mean(as.numeric(sampleratio), na.rm=T),
                     meanMSE=mean(MSE, na.rm=T),
                     LCIMSE=meanMSE-miscset:::confint.numeric(MSE, na.rm = T),
                     UCIMSE=meanMSE+miscset:::confint.numeric(MSE, na.rm = T),
                     file=index)
  return(res_sum)
}

# For normal running
{
  #for dataset=1
  performance_1=lapply(1:1, function(seed_run) {
    miss=30; corr=0.52
    incomplete_feature = list(miss=miss, corr=corr, datacode= 1, seed = seed_run, GA=T, cutratio=2, datatype="incomplete")
    # print(unlist(incomplete_feature))
    incomplete_res = mfit(incomplete_feature)
    #if(getres >0 ){
    both_feature= list(miss=miss, corr=corr, datacode= 1, seed = seed_run, GA=T, cutratio=2, datatype="both")
    both_res=mfit(both_feature)
    #print(both_res)
    #both_res
    # #}else{both_res=NULL}
    finalres=rbind(incomplete_res,both_res)
    #print(finalres)
    finalres
  })
  #print(performance_1)
  performance_1_df=do.call(rbind, performance_1)
  summarizer_bay(performance_1_df)
  write.csv(performance_1_df, "performance_1_df.csv")

  #for dataset=2
  performance_2=lapply(1:1, function(seed_run) {
    miss=10; corr=0.52
    incomplete_feature= list(miss=miss, corr=corr, datacode= 2, seed = seed_run, GA=T, cutratio=2, datatype="incomplete")
    # print(unlist(incomplete_feature))
    incomplete_res=mfit(incomplete_feature)
    #if(getres >0 ){
    both_feature= list(miss=miss, corr=corr, datacode= 2, seed = seed_run, GA=T, cutratio=2, datatype="both")
    both_res=mfit(both_feature)
    #print(both_res)
    #both_res
    # #}else{both_res=NULL}
    finalres=rbind(incomplete_res,both_res)
    #print(finalres)
    finalres
  })
  performance_2_df=do.call(rbind, performance_2)
  summarizer_bay(performance_2_df)
  write.csv(performance_2_df, "performance_2_df.csv")
  
  #for dataset=3
  performance_3=lapply(1:1, function(seed_run) {
    miss=30; corr=0.52
    incomplete_feature= list(miss=miss, corr=0.52, datacode= 4, seed = seed_run, GA=T, cutratio=2, datatype="incomplete")
    # print(unlist(incomplete_feature))
    incomplete_res=mfit(incomplete_feature)
    # if(getres >0 ){
    both_feature= list(miss=miss, corr=0.52, datacode= 4, seed = seed_run, GA=T, cutratio=2, datatype="both")
    both_res=mfit(both_feature)
    #print(both_res)
    #both_res
    # #}else{both_res=NULL}
    finalres=rbind(incomplete_res,both_res)
    #print(finalres)
    finalres
  })
  performance_3_df=do.call(rbind, performance_3)
  summarizer_bay(performance_3_df)
  write.csv(performance_3_df, "performance_3_df.csv")
  
  #for dataset=4
  performance_4=lapply(1:1, function(seed_run) {
    incomplete_feature= list(miss=10, corr=0.52, datacode= 5, seed = seed_run, GA=T, cutratio=2, datatype="incomplete")
    # print(unlist(incomplete_feature))
    incomplete_res=mfit(incomplete_feature)
    #if(getres >0 ){
    both_feature= list(miss=10, corr=0.52, datacode= 5, seed = seed_run, GA=T, cutratio=2, datatype="both")
    both_res=mfit(both_feature)
    #print(both_res)
    #both_res
    # #}else{both_res=NULL}
    finalres=rbind(incomplete_res,both_res)
    print(finalres)
    finalres
  })
  performance_4_df=do.call(rbind, performance_4)
  write.csv(performance_4_df, "performance_4_df.csv")
  
  
  for(seed_run in 1:30){
    feature= list(50,0.8, seed_run, 63, GA=T, cutratio=5)
    performance=fit(feature)
  }
}
