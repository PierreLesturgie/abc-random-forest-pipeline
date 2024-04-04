################################################################
#----------------  PIPELINE RANDOM FORESTS ABC ----------------#
################################################################


## ****** author  = Pierre Lesturgie ****** ##
## ****** version = 0.2 ****** ##
## ****** email = pierrelesturgie@outlook.fr ****** ##
## ****** date = 04.04.2024 ****** ##


#------- MODEL SELECTION ------#

# >>>>>> MS-1: model.selection.random.forest() 
# performs  model selection from a target + a named list of summary statistics simulated from different scenarios

# ==> Input : 
#       *** target : observed summary statistics
#       *** list_sumstat : Named list of summary statistics computed from simulated models
#       *** dir : directory for output files
#       *** analysis : "LDA", "basic" or "all", to compute the ms adding linear discriminant analysis axes as sumstat, without, or both. 
#       *** ntree : number of trees used to grow the forest
#       *** predictions : does the function computes model choice (or only grows forest) ?
#       *** compute_oob : does the function computes the evolution of out-of-bag error with the number of trees ? 
#       *** importance_variable : does the function tests the importance of the variables in the decision ? 
#       *** subset : number of observation to subset each dataframe. 
#       *** LDA_plot : does the function computes the LDA plot (very similar to a pca)
# ==> Output : 
#       *** returns a list with the confusion matrix and the output of analyses asked to the function. 


model.selection.random.forest<-function(target, list_sumstat,directory=c(""),analysis=c("all"),ntree=500,predictions=T,
                                        compute_oob=T,importance_variable=T,subset=NULL,LDA_plot=T,
                                        paral = F){
  require("abcrf")
  
  prep = prepa.sumstat(list_sumstat)
  
  if (length(target)!=ncol(prep$sumstat)){stop(paste0("Not the same number of statistics in target and simulated data!"))}
  
  y = as.factor(prep$y)
  sumstat = prep$sumstat
  data_rf<-data.frame(y,sumstat)
  colnames(target)<-colnames(data_rf)[-1]
  
  if(analysis=="LDA"){
    cat("Building random forest with Linear Discriminant Analysis...","\n")
    rfLDA<-abcrf(y~., data = data_rf, lda=T, ntree=ntree)
    cat("...Done.","\n")
    print(rfLDA)
    cat("\n","\n")
    if (predictions==T){
      cat("Model prediction with Linear Discriminant Analysis...","\n")
      pLDA<-predict(rfLDA, target, training = data_rf, paral = paral,ntree = ntree)
      cat("...Done.","\n")
      print(pLDA)
      ppLDA<-cbind(pLDA$vote,pLDA$post.prob);colnames(ppLDA)<-c(colnames(pLDA$vote),"post.prob")
      write.table(ppLDA,paste0(directory,"/predictions_LDA.txt"))
      cat("\n","\n")
    }
    if(compute_oob==T){
      cat("Out-of-bag error analysis...","\n")
      oob_rfLDA<-err.abcrf(rfLDA,data_rf)
      pdf(paste0(directory,"/OOB_LDA.pdf"))
      plot(oob_rfLDA)
      dev.off()
      cat("...Done.","\n")
      cat("\n")
    }
    if(importance_variable==T){
      impLDA<-variableImpPlot(object = rfLDA)
    }
  }else if (analysis=="all"){
    cat("Building random forest...","\n")
    rf<-abcrf(y~., data = data_rf, lda=F, ntree=ntree)
    cat("...Done.","\n")
    print(rf)
    if (predictions==T){
      cat("Model prediction...","\n")
      p<-predict(rf, target, training = data_rf, paral = paral,ntree = ntree)
      cat("...Done.","\n")
      print(p)
      pp<-cbind(p$vote,p$post.prob);colnames(pp)<-c(colnames(p$vote),"post.prob")
      write.table(pp,paste0(directory,"/predictions.txt"))
      cat("\n","\n")
    }
    if(compute_oob==T){
      cat("Out-of-bag error analysis...","\n")
      oob_rf<-err.abcrf(rf,data_rf)
      pdf(paste0(directory,"/OOB.pdf"))
      plot(oob_rf)
      dev.off()
      cat("...Done.","\n")
      cat("\n")
    }
    if(importance_variable==T){
      imp<-variableImpPlot(object = rf)
    }
    cat("Building random forest with Linear Discriminant Analysis...","\n")
    rfLDA<-abcrf(y~., data = data_rf, lda=T, ntree=ntree)
    cat("...Done.","\n")
    print(rfLDA)
    cat("\n","\n")
    if (predictions==T){
      cat("Model prediction with Linear Discriminant Analysis...","\n")
      pLDA<-predict(rfLDA, target, training = data_rf, paral = paral,ntree = ntree)
      cat("...Done.","\n")
      print(pLDA)
      ppLDA<-cbind(pLDA$vote,pLDA$post.prob);colnames(ppLDA)<-c(colnames(pLDA$vote),"post.prob")
      write.table(ppLDA,paste0(directory,"/predictions_LDA.txt"))
      cat("\n","\n")
    }
    if(compute_oob==T){
      cat("Out-of-bag error analysis...","\n")
      oob_rfLDA<-err.abcrf(rfLDA,data_rf)
      pdf(paste0(directory,"/OOB_LDA.pdf"))
      plot(oob_rfLDA)
      dev.off()
      cat("...Done.","\n")
      cat("\n")
    }
    if(importance_variable==T){
      impLDA<-variableImpPlot(object = rfLDA)
    }
  }else if (analysis == "basic"){
    cat("Building random forest...","\n")
    rf<-abcrf(y~., data = data_rf, lda=F, ntree=ntree)
    cat("...Done.","\n")
    print(rf)
    cat("\n")
    if(predictions==T){
      cat("Model prediction...","\n")
      p<-predict(rf, target, training = data_rf, paral = paral,ntree = ntree)
      cat("...Done.","\n")
      print(p)
      pp<-cbind(p$vote,p$post.prob);colnames(pp)<-c(colnames(p$vote),"post.prob")
      write.table(pp,paste0(directory,"predictions.txt"))
      cat("\n","\n")
    }
    if(compute_oob==T){
      cat("Out-of-bag error analysis...","\n")
      oob_rf<-err.abcrf(rf,data_rf)
      cat("...Done.","\n","\n")
    }
    if(importance_variable==T){
      imp<-variableImpPlot(object = rf)
    }
  }
  
  
  if (LDA_plot==T){
    y2<-as.factor(c(prep$y,"TARGET"))
    sumstat2<-rbind(sumstat,target)
    data_rf2<-data.frame(y2,sumstat2)
    rf2<-abcrf(y2~., data = data_rf2, 
               lda=T, ntree=ntree,
               paral=paral)
    pdf(paste0(directory,"/LDA_plot.pdf"))
    plot(rf2,training=data_rf2)
    
    dev.off()
  }
  
  
  i=1
  if (analysis=="all"){
    final<-list()
    finalLDA<-list()
    final[[i]]<-rf$model.rf$confusion.matrix
    final[[i+1]]<-rf$model.rf$prediction.error
    finalLDA[[i]]<-rfLDA$model.rf$confusion.matrix
    finalLDA[[i+1]]<-rfLDA$model.rf$prediction.error
    names(final)[1:2]<-c("confusion_matrix","Prior_error_rate")
    names(finalLDA)[1:2]<-c("confusion_matrix","Prior_error_rate")
    i<-i+2
    if(predictions==T){
      final[[i]]<-p
      finalLDA[[i]]<-pLDA
      names(final)[i]<-c("predictions")
      names(finalLDA)[i]<-c("predictions")
      i=i+1
    }
    if(compute_oob==T){
      final[[i]]<-oob_rf
      finalLDA[[i]]<-oob_rfLDA
      names(final)[i]<-c("OOB_ntrees")
      names(finalLDA)[i]<-c("OOB_ntrees")
      i=i+1
    }
    if(importance_variable==T){
      final[[i]]<-imp
      finalLDA[[i]]<-impLDA
      names(final)[i]<-c("Variable_importance")
      names(finalLDA)[i]<-c("Variable_importance")
      i=i+1
    }
    ALL_FINAL<-list(final,finalLDA);names(ALL_FINAL)<-c("No_LDA","LDA")
    return(ALL_FINAL)
  }else if(analysis=="LDA"){
    finalLDA<-list()
    finalLDA[[i]]<-rfLDA$model.rf$confusion.matrix
    finalLDA[[i+1]]<-rfLDA$model.rf$prediction.error
    names(finalLDA)[1:2]<-c("confusion_matrix","Prior_error_rate")
    i<-i+2
    if(predictions==T){
      finalLDA[[i]]<-pLDA
      names(finalLDA)[i]<-c("predictions")
      i=i+1
    }
    if(compute_oob==T){
      finalLDA[[i]]<-oob_rfLDA
      names(finalLDA)[i]<-c("OOB_ntrees")
      i=i+1
    }
    if(importance_variable==T){
      finalLDA[[i]]<-impLDA
      names(finalLDA)[i]<-c("Variable_importance")
      i=i+1
    }
    return(finalLDA)
  }else if(analysis=="basic"){
    final<-list()
    final[[i]]<-rfLDA$model.rf$confusion.matrix
    final[[i+1]]<-rfLDA$model.rf$prediction.error
    names(final)[1:2]<-c("confusion_matrix","Prior_error_rate")
    i<-i+2
    if(predictions==T){
      final[[i]]<-p
      names(final)[i]<-c("predictions")
      i=i+1
    }
    if(compute_oob==T){
      final[[i]]<-oob_rf
      names(final)[i]<-c("OOB_ntrees")
      i=i+1
    }
    if(importance_variable==T){
      final[[i]]<-imp
      names(final)[i]<-c("Variable_importance")
      i=i+1
    }
    return(final)
  }
}

# >>>>>> MS-2: pcabc()
# performs a PCA with simulated and target summary statistics 

# ==> Input : 
#       *** models : names of models to compute (only important for loading sumstats and for output files)
#       *** dir : directory for INPUT and output files
#           --> sumstat files named "sumstat'model'.txt" (ex sumstat_FIM.txt) << only the folded-SFS >>
#           --> target file named "target.txt" << only the folded-SFS >>
#       *** nind : number of individuals
#       *** subset : number of observation to subset each dataframe. 
#       *** pcs : number of pcs (if NULL will be asked) 
# ==> Output : 
#       *** write pdf of the pca
#       *** returns a pca object

pcabc<-function(target,list_sumstat,directory=c("./"),subset=NULL,pcs=NULL){
  require(ade4)
  
  prep = prepa.sumstat(list_sumstat,subset=subset)
  
  if (length(target)!=ncol(prep$sumstat)){stop(paste0("Not the same number of statistics in target and simulated data!"))}
  
  y<-as.vector(c(prep$y,"TARGET"))
  colnames(target)<-colnames(prep$sumstat)
  
  sumstat<-rbind(prep$sumstat,target)
  
  data_rf<-data.frame(y,sumstat)
  
  if (is.null(pcs)==F){
    pca<-dudi.pca(df = data_rf[, -1], scannf = FALSE, nf = pcs)
  }
  pca<-dudi.pca(data_rf[,-1])
  
  col<-c()
  for (i in 1:length(list_sumstat)){
    col=c(col,rep(i,prep$nsims[i]))
  }
  
  pcs<-pca$nf
  d<-(pca$eig/sum(pca$eig))*100
  
  pdf(paste0(directory,"pca.pdf"),onefile = T)
  for (j in 1:(pcs-1)){
    plot(pca$tab[-1,c(1,j+1)],col=col,pch = 19,cex=0.2,xlab=paste0("PC1 (",round(d[1]),"%)"),
         ylab=paste0("PC",j+1," (",round(d[j+1]),"%)")) + points(pca$tab[1,c(1,j+1)],
                                                                 col=i+1,pch = 18,cex=1)
    legend("topleft",legend = c(names(list_sumstat),"Target"),col=c(1:(i+1)),pch = 19)
    
  }
  
  dev.off()
  return(pca)
}

# >>>>>> MS-3-hide: prepa.sumstat()
# simply builds dataset for model selection. Returns a list with sumstat df + vector of models

prepa.sumstat=function(list_sumstat,subset=NULL){
  if(is.null(names(list_sumstat))){stop(paste0("Need to name the list 'list_sumstat' by the different models tested"))}
  
  sumstat<-c()
  nsims<-c()
  for (i in 1:length(list_sumstat)){
    temp<-list_sumstat[[i]]
    if(is.null(subset)==F){
      m<-c(sample(x = c(1:nrow(temp)),size = subset,replace = F))
      temp<-temp[m,]
    }
    sumstat<-rbind(sumstat,temp)
    nsims<-rbind(nsims,nrow(temp))
  }
  
  y<-c()
  for (i in 1:length(list_sumstat)){
    y<-rbind(y,as.matrix(rep(names(list_sumstat)[i],nsims[i])))
  }
  y<-as.vector(y)
  res=list(sumstat,y,nsims); names(res) = c("sumstat","y","nsims")
  return(res)
}


#------- PARAMETER ESTIMATION -------#

# >>>>>> PE-1 : param.estimation.random.forest()

# Performs parameter estimation, and optionally cross-validation

# ==> Input : 
#       *** sumstat : data frame w/ sumstat data (SFS + other wanted sumstats)
#       *** priors : priors corresponding to the sumstat data
#       *** target : targetted sumstat (need to have exactly the same stats computed than sumstat)
#       *** intv_qtil : steps between quantiles (for the plots)
#       *** nval : number of points used to plot the values from the quantiles
#       *** n_tree : number of trees used to grow the forest
#       *** CI : bounds of the confidence interval
#       *** dir : directory for output files
#       *** cross_validation : TRUE for computation of abcRF package implemented cross validation <<< very time consuming >>>
#       *** densityPlot : TRUE to compute and save the abcRF implemented density plot
#       *** save_rf : TRUE to save the random forests (only if cross_validation=F) <<< very heavy objects >>>
# ==> Output : 
#       *** print graphs in dir
#       *** print tables in dir
#       *** returns a list with
#           --> a list of reference tables built for the random forests
#           --> a list of estimated values data frames
#           --> a list of out-of-bag error values according to the number of trees for each forest
#           --> the formatted dataframe used to plot
#           --> a list of prediction of the out-of-bag error for each forest (cross-val, only if cross_validation=T)
#           --> a list of the random forest objects (only if save_rf=T)



param.estimation.random.forest<-function(sumstat,priors,target,n_tree=500,dir="",cross_validation=T,predict_OOB_CV=F,
                                         intv_qtil=1e-4,nval=1000,npods=1000,CI=c(0.025,0.975),paral=F,
                                         densityPlot=T,save_rf=F,write_cv=F){
  
  require(abcrf)
  
  q<-seq(0, 1, by = intv_qtil)
  q<-as.vector(q)
  cross_val_output<-c()
  cross_val<-list()
  cross_val_result = c()
  error<-list()
  data_list<-list()
  prediction_list<-list()
  stat_table<-c()
  
  if (predict_OOB_CV==F){
    if (save_rf==T){rf_list<-list()}
  }
  
  for (i in 1:ncol(priors)){
    cat("==> Analysis of variable",colnames(priors)[i],"\n")
    data_list[[i]]<-as.data.frame(cbind(priors[,i],sumstat))
    colnames(data_list[[i]])[1]<-"Temp"
    cat("     >> building the random forest...","\n")
    rf<-regAbcrf(Temp~.,data_list[[i]], ntree = n_tree, paral = paral)
    cat("     Done.","\n")
    cat("********* Info random forest on parameter",colnames(priors)[i],":","\n")
    print(rf)
    cat("\n","\n")
    cat("     >> computing estimation...","\n")
    prediction_list[[i]]<-predict(rf, target,data_list[[i]], 
                                  paral = paral, quantiles = q)
    cat("     Done.","\n","\n")
    stat_temp<-cbind(prediction_list[[i]]$expectation, prediction_list[[i]]$med, prediction_list[[i]]$variance, 
                     prediction_list[[i]]$variance.cdf, prediction_list[[i]]$post.NMAE.mean)
    colnames(stat_temp)<-c("mean", "median", "variance", "variance.cdf", "post.NMAE.mean")
    stat_table<-rbind(stat_table,stat_temp)
    
    cat("     >> computing out-of-bag MSE...","\n")
    error[[i]]<-err.regAbcrf(rf,training = data_list[[i]],paral=paral)
    pdf(paste0(dir,"OOB_",colnames(priors)[i],".pdf"))
    plot(error[[i]])
    dev.off()
    cat("     Done.","\n","\n")
    
    if (predict_OOB_CV==T){
      cat("     >> computing cross-validation...","\n")
      cross_val[[i]]<-predictOOB(rf,training = data_list[[i]],paral=paral)
      temp_cv<-data.frame(cross_val[[i]]$MSE,cross_val[[i]]$NMAE,
                          cross_val[[i]]$MSE.med,cross_val[[i]]$NMAE.med,cross_val[[i]]$coverage)
      colnames(temp_cv)<-c("OOB-MS-mean","OOB-NMAE-mean","OOB-MS-median","OOB-NMAE-median","CI-coverage")
      cross_val_output<-rbind(cross_val_output,temp_cv)
      cat("     Done.","\n","\n")
      
    }
    
    if (cross_validation==T){
      cat("     >> computing cross-validation...","\n")
      cv = cross.validation.rf(rf=rf,ref_table=data_list[[i]],priors=priors[,i],var=colnames(priors)[i],
                          sumstats=sumstat,quantiles = q,intv_qtil=intv_qtil,nval=nval,npods=npods,CI=c(0.025,0.975),write_cv = F,
                          paral=paral)
      cat("     Done.","\n","\n")
      cross_val_result = cbind(cross_val_result,cv)
    }
    
    if (densityPlot==T){
      cat("     >> printing abcrf density plot...","\n")
      pdf(paste0(dir,"densityPlot_abcrf_",colnames(priors)[i],".pdf"))
      densityPlot(rf,obs = target,training = data_list[[i]],paral=paral)
      dev.off()
      cat("     Done.","\n","\n")
    }
    
    if (predict_OOB_CV==F){
      if (save_rf==T){
        rf_list[[i]]<-rf;names(rf_list)[i]<-colnames(priors)[i]
      } 
    }
    rm(rf)
  }
  
  if (predict_OOB_CV==T){
    names(cross_val)=colnames(priors)
    write.table(cross_val_output,"cross_val.txt")
  }
  
  names(prediction_list)=names(data_list)=names(error)=colnames(priors)
  
  df<-formatting.RF.est(prediction_list,dir = dir,byqtil=intv_qtil,nval=nval,CI=NULL)
  
  estimations<-estimate.plot.RF(df,priors,save = T,dir = dir,titles_param_plot = NULL)
  
  # For confidence interval
  est<-param.estimates(df=df,CI=CI)
  
  write.table(est,file = paste0(dir,"estimates.txt"))
  write.table(df,file = paste0(dir,"reg.txt"))
  if (cross_validation==T){write.table(cross_val_result,paste0(dir,"cross_val.txt"))}
  
  
  
  #e<-err.regAbcrf(r.rf_fim_cuv_c, Nm_fim_cuv_c, paral = T)
  
  if (save_rf==T){
    fin<-list(data_list,rf_list,est,error,df,cross_val_result,cross_val);names(fin)<-c("ref_table","random_forests","estimated_values","error","reg_df")
  }else{
    fin<-list(data_list,est,error,df,cross_val_result,cross_val);names(fin)<-c("ref_table","estimated_values","error","reg_df","cross_validation","predict_OOB_CV")
  }
  
  return(fin)
  
}





# >>>>>> PE-2-hide : formatting.RF.est()
# setting up data base for parameters estimation 

# ==> Input : 
#       *** data_prediction : list of quantiles of parameter estimations
#       *** byqtil : interval within quantiles chosen for giving q quantiles in regAbcrf() function
#       *** nval : number of quantiles/parameters (for graphic purposes, "resolution")
#       *** dir : directory for output files
#       *** CI : bounds of the confidence interval
# ==> Output : 
#       *** write table of CI for parameters
#       *** df_fin: a data frame of nval rows and ncol variables corresponding to estimated values for each variables. 


formatting.RF.est<-function(data_prediction, dir,byqtil=1e-4,nval=1000,CI=c(0.025,0.975)){
  q<-seq(from=0, to=1, by=byqtil)
  b<-seq(from=1, to=1/byqtil, by=1/(byqtil*nval))
  
  df_temp<-c()
  for (i in 1:length(data_prediction)){
    data_temp<-as.data.frame(data_prediction[[i]]$quantiles)
    colnames(data_temp)<-c(q);data_temp<-t(data_temp)
    df_temp<-cbind(df_temp,data_temp)
  }
  df_temp<-cbind(df_temp,q)
  df_fin<-df_temp[b,]
  
  rownames(df_fin)<-NULL
  colnames(df_fin)<-c(names(data_prediction),"Quantiles")
  df_fin<-as.data.frame(df_fin)
  
  
  if (is.null(CI)==F){
    # For confidence intervalle
    CInames<-paste0(as.character(CI*100),"%")
    CI<-c(which(df_fin$Quantiles>=CI[1])[1],which(df_fin$Quantiles>=CI[2])[1])
    if(is.na(CI[2])){CI[2]<-nrow(df_fin)}
    ic<-df_fin[CI,]
    ic<-t(ic)
    colnames(ic)<-c("2.5%","97.5%")
    ic<-ic[-nrow(ic),]
    
    mediane<-c()
    mode<-c()
    for (i in 1:length(data_prediction)){
      mediane<-c(mediane,data_prediction[[i]]$med)
      temp<-density(df_fin[,i]);mode_temp<-temp$x[which(temp$y==max(temp$y))]
      mode<-c(mode,mode_temp)
    }
    est<-data.frame(mode,mediane,ic[,1],ic[,2]);colnames(est)<-c("Mode","Median",CInames[1],CInames[2])
    
    write.table(file = paste0(dir,"estimates.txt"),est)
  }
  
  
  return(df_fin)
}

# >>>>>> PE-3-hide : estimate.plot.RF()
# parameters estimation plot

#       <<< tranformation of quantiles weights to find back approximation of the real posterior distribution : >>>
#       <<< quantile[i]/(param[i]-param[i-1]) ==> gives a 2*(nval-1)*nparam output list >>>

# ==> Input : 
#       *** data : data frame w/ quantiles corressponding (output from formatting.RF.est())
#       *** priors : priors corresponding to the data.frame
#       *** dir : directory for output files
#       *** titles_param_plot : vector names for the plots (dim is the number of variables)
#       *** x_lim : TRUE to center the plot on the max values of prior (false centered on estimates)
# ==> Output : 
#       *** print graphs in dir
#       *** returns a list of the dataframes to each variable plotted + priors. 

estimate.plot.RF<-function(data, priors, save=T, dir="./",titles_param_plot=NULL,x_lim=F){
  
  param_estim<-list()
  for (j in 1:(ncol(data)-1)){
    d<-c()
    for (i in 2:(length(data[,1]))){
      temp<-c(data[i,j]-data[i-1,j])
      d<-c(d,temp)
    }
    if (length(which(d==0))!=0){
      d[which(d==0,arr.ind = T)]<-min(d[-(which(d==0,arr.ind = T))])
    }
    
    d<-c(d,data[1,j])
    temp2<-data[,j]
    q<-data[,ncol(data)]/d
    param_estim[[j]]<-cbind(temp2,q);colnames(param_estim[[j]])<-c(colnames(data)[j],"Quantiles")
  }
  names(param_estim)<-colnames(data)[1:(ncol(data)-1)]
  
  my.plot.here<-function(datum,prior,titles_param_plot=NULL,x_lim=T){
    
    if (is.null(titles_param_plot)==F){title<-titles_param_plot}else {tit=''}
    if (x_lim==T){
      p<-plot(density(datum[,1]), main = tit, xlab=colnames(datum)[1],xlim=c(min(prior),max(prior))) + 
        lines(density(prior), col = "grey") + 
        lines(density(datum[,1]))
    }else{
      p<-plot(density(datum[,1]), main = tit, xlab=colnames(datum)[1]) + 
        lines(density(prior), col = "grey") + 
        lines(density(datum[,1]))
    }
    #print(p)
  }
  
  if (save==T){
    pdf(paste0(dir,"estimation.pdf"),onefile = T)
    for (i in 1:(ncol(data)-1)){
      my.plot.here(datum = param_estim[[i]],prior = priors[,i],x_lim = x_lim)
    }
    dev.off()
  }else{
    for (i in 1:(ncol(data)-1)){
      my.plot.here(datum = param_estim[[i]],prior = priors[,i],x_lim = x_lim)
    }
  }
  
  
  
  fin<-param_estim ; fin[[length(param_estim)+1]]<-priors
  names(fin)<-c(names(param_estim),"priors")
  return(fin)
}

# >>>>>> PE-4-hide : param.estimates()
# extract parameter estimates and CI

# ==> Input : 
#       *** df : data frame of quantiles for all estimates
#       *** CI : bounds of the confidence interval
# ==> Output : 
#       *** data frame of estimates (mode, median, CI)

param.estimates<-function(df,CI=c(0.025,0.975)){
  CInames<-paste0(as.character(CI*100),"%")
  CI<-c(which(df$Quantiles>=CI[1])[1],which(df$Quantiles>=CI[2])[1])
  if(is.na(CI[2])){CI[2]<-nrow(df)}
  ic<-df[CI,-ncol(df)]
  if (class(ic)!="data.frame"){ic<-as.data.frame(ic)}
  rownames(ic)<-CInames;colnames(ic)<-colnames(df)[-ncol(df)]
  
  mediane<-df[which(df$Quantiles==0.5),-ncol(df)]
  if (class(mediane)!="data.frame"){mediane<-as.data.frame(mediane)}
  rownames(mediane)<-"median";colnames(mediane)<-colnames(df)[-ncol(df)]
  mode<-c()
  for (i in 1:(ncol(df)-1)){
    temp<-density(df[,i]);mode_temp<-temp$x[which(temp$y==max(temp$y))]
    mode<-c(mode,mode_temp)
  }
  mode<-t(as.data.frame(mode));colnames(mode)<-colnames(mediane)
  est<-rbind(mode,mediane,ic)
  return(est)
}

# >>>>>> PE-5-hide : cross.validation.rf()
# performs a *custom* cross-validation 
# procedure is simply to perform prediction based on Pseudo-Observed Dataset (POD)

cross.validation.rf<-function(rf, ref_table,npods=1000,priors,sumstats,var,
                              quantiles,intv_qtil=intv_qtil,CI=c(0.025,0.975),paral=F,write_cv=F,nval=1000){
  
  priors<-data.frame(priors)
  stats = pods = c()
  pb <- txtProgressBar(min = 0, max = npods, style = 3)
  for(i in 1:npods){
    idx<-sample(nrow(priors),size = 1)
    pod<-priors[idx,]
    target<-sumstats[idx,]
    #print(target)
    target<-data.frame(target)
    #colnames(target)<-colnames(ref_table)[-1]
    
    pred<-predict(object=rf, obs=target,training=ref_table,paral = paral,quantiles = quantiles)
    list_pred<-list(pred);names(list_pred)<-"temp"
    df<-formatting.RF.est(list_pred,byqtil=intv_qtil,nval=nval,CI=NULL)
    estimations<-estimate.plot.RF(df,priors,save = FALSE,titles_param_plot = NULL)
    
    est<-param.estimates(df=df,CI=CI)
    est<-est[-2,]; est<-t(est)
    #print(est)
    stat_temp<-cbind(pred$expectation, pred$med, est[1,1], pred$variance, 
                     pred$variance.cdf, pred$post.NMAE.mean,est[1,2],est[1,3])
    
    
    colnames(stat_temp)<-c("mean", "median", "mode","variance", "variance.cdf", "post.NMAE.mean","quantile=0.025", "quantile=0.975")
    stats = rbind(stats,stat_temp)
    pods = rbind(pods,pod)
    
    
    if(write_cv==T){
      write.table(stat_temp,file=paste0(dir,"/",var,"_cv_stats_",i,".txt"))
      write.table(pod,file=paste0(dir,"/",var,"_cv_pod_",i,".txt"))
      write.table(target,file=paste0(dir,"/",var,"_cv_ss_",i,".txt"))
    }
    
    
    setTxtProgressBar(pb, i)
    
  }
  pa = post.analysis.cv(param = pods,stats = stats)
  colnames(pa) = var
  return(pa)
}


# >>>>>> PE-6-hide : post.analysis.cv()
# performs coverage SME and SRMSE from PODs and estimates in cross val
post.analysis.cv = function(param,stats){
  #### STATS : "mean" "median" "mode" "variance" "variance.cdf" "post.NMAE.mean" "quantile=0.025" "quantile=0.975"
  
  # coverage = 95% coverage (CI) Number of observation included in the CI
  # SME = Squared Mean Error, as in Bruno A. Walther and Joslin L. Moore, Ecography 2005
  # SRMSE = scaled root mean square error, as in Bruno A. Walther and Joslin L. Moore, Ecography 2005
  
  sme = function(obs,pods){return(sum((pods-obs)/ obs)/length(obs))}
  srmse = function(obs,pods){return((sum(((pods-obs)/obs)^2)/length(obs))^0.5)} 
  coverage = function(obs,pod_low_ci,pod_high_ci){
    coverage<-0
    for (j in 1:length(obs)) {
      if ((obs[j])<= pod_high_ci[j] && (obs[j])>=pod_low_ci[j]){
        coverage<-coverage+1
      }  
    }
    coverage<-coverage/length(obs)
    return(coverage)
  }
  
  cov = coverage(obs = param[,1],pod_low_ci = stats[,7],pod_high_ci = stats[,8])
  mean_SME = sme(obs = param[,1], pods = stats[,1])
  mean_SRMSE = srmse(obs = param[,1], pods = stats[,1])
  mediane_SME = sme(obs = param[,1], pods = stats[,2])
  mediane_SRMSE = srmse(obs = param[,1], pods = stats[,2])
  mode_SME = sme(obs = param[,1], pods = stats[,3])
  mode_SRMSE = srmse(obs = param[,1], pods = stats[,3])
  
  finale=rbind(cov,mean_SME,mean_SRMSE,mediane_SME,mediane_SRMSE,mode_SME,mode_SRMSE)
  rownames(finale)<-c("coverage95%","mean_SME","mean_SRMSE","mediane_SME","mediane_SRMSE","mode_SME","mode_SRMSE")
  return(finale)
}





