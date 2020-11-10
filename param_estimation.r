################################################################
#------ FUNCTIONS PARAM ESTIMATION W/ RANDOM FORESTS ABC ------#
################################################################


#-------function 1 ==> Formatting data in order to compute the abcRF regression : formatting.RF.est()



# >>>>>> function 1 : data base form for parameters estimation 

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
    est<-data.frame(ic[,1],mode,mediane,ic[,2]);colnames(est)<-c(CInames[1],"Mode","Median",CInames[2])
  
    write.table(file = paste0(dir,"/estimates.txt"),est)
  }
  
  
  return(df_fin)
}



# >>>>>> function 2 : parameters estimation plot####

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

estimate.plot.RF<-function(data, priors, save=T, dir="",titles_param_plot=NULL,x_lim=F){
  
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
    q<-data[,4]/d
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
    print(p)
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



# >>>>>> function 3 : do the whole estimation from sumstat data 

#       <<< Needs formatting.RF.est() & estimate.plot.RF() >>>

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


param.estimation.random.forest<-function(sumstat,priors,target,intv_qtil=1e-4,nval=1000,n_tree=500,CI=c(0.025,0.975),dir="",cross_validation=T,paral=T,densityPlot=T,save_rf=F){
  
  require(abcrf)
  
  q<-seq(0, 1, by = intv_qtil)
  q<-as.vector(q)
  
  cross_val<-list()
  error<-list()
  data_list<-list()
  prediction_list<-list()
  stat_table<-c()
  if (cross_validation==F){
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
    pdf(paste0(dir,"/OOB_",colnames(priors)[i],".pdf"))
    plot(error[[i]])
    dev.off()
    cat("     Done.","\n","\n")
    
    if (cross_validation==T){
      cat("     >> computing cross-validation...","\n")
      cross_val[[i]]<-predictOOB(rf,training = data_list[[i]],paral=paral)
      cat("     Done.","\n","\n")
      
      
    }
    if (densityPlot==T){
      cat("     >> printing abcrf density plot...","\n")
      pdf(paste0(dir,"/densityPlot_abcrf_",colnames(priors)[i],".pdf"))
      densityPlot(rf,obs = target,training = data_list[[i]],paral=paral)
      dev.off()
      cat("     Done.","\n","\n")
    }
    if (cross_validation==F){
      if (save_rf==T){
        rf_list[[i]]<-rf;names(rf_list)[i]<-colnames(priors)[i]
      } 
    }
    rm(rf)
  }
  if (cross_validation==T){names(cross_val)=colnames(priors)}
  
  names(prediction_list)=names(data_list)=names(error)=colnames(priors)
  
  df<-formatting.RF.est(prediction_list,dir = dir,byqtil=intv_qtil,nval=nval,CI=NULL)
  
  estimations<-estimate.plot.RF(df,priors,save = T,dir = dir,titles_param_plot = NULL)
  
  # For confidence intervalle
  CInames<-paste0(as.character(CI*100),"%")
  CI<-c(which(df$Quantiles>=CI[1])[1],which(df$Quantiles>=CI[2])[1])
  if(is.na(CI[2])){CI[2]<-nrow(df)}
  ic<-df[CI,]
  ic<-t(ic)
  colnames(ic)<-c("2.5%","97.5%")
  ic<-ic[-nrow(ic),]
  
  mediane<-c()
  mode<-c()
  for (i in 1:length(prediction_list)){
    mediane<-c(mediane,prediction_list[[i]]$med)
    temp<-density(df[,i]);mode_temp<-temp$x[which(temp$y==max(temp$y))]
    mode<-c(mode,mode_temp)
  }
  
  est<-data.frame(ic[,1],mode,mediane,ic[,2]);colnames(est)<-c(CInames[1],"Mode","Median",CInames[2])
  
  write.table(file = paste0(dir,"/estimates.txt"),est)
  
  #e<-err.regAbcrf(r.rf_fim_cuv_c, Nm_fim_cuv_c, paral = T)
  
  if (cross_validation==F) {
    if (save_rf==T){
      fin<-list(data_list,rf_list,est,error,df);names(fin)<-c("ref_table","random_forests","estimated_values","error","reg_df")
    }else{
      fin<-list(data_list,est,error,df);names(fin)<-c("ref_table","estimated_values","error","reg_df")
    }
    
  } else{fin<-list(data_list,est,error,cross_val,df);names(fin)<-c("ref_table","estimated_values","error","cross_validation","reg_df")}
  
  
  return(fin)
  
}

# >>>>>> function 3bis : do the whole estimation from RAW sumstat data (useful for the cluster)

#       <<< Needs formatting.RF.est() & estimate.plot.RF() >>>
#       <<< Needs functions from Stefano to compute mpd and TD >>>
#       <<< Needs a folder with sumstats, priors and target >>>

# ==> Input : 
#       *** dir : directory for INPUT and output files
#           --> sumstat file named "sumstat.txt" << only the folded-SFS >>
#           --> priors file named "priors.txt"
#           --> target file named "target.txt" << only the folded-SFS >>
#       *** Nindiv : number of diploid individuals
#       *** intv_qtil : steps between quantiles (for the plots)
#       *** nval : number of points used to plot the values from the quantiles
#       *** n_tree : number of trees used to grow the forest
#       *** CI : bounds of the confidence interval
#       *** cross_validation : TRUE for computation of abcRF package implemented cross validation <<< very time consuming >>>
#       *** densityPlot : TRUE to compute and save the abcRF implemented density plot
# ==> Output : 
#       *** print graphs in dir
#       *** print tables in dir
#       *** returns a list with
#           --> a list of reference tables built for the random forests
#           --> a list of estimated values data frames
#           --> a list of out-of-bag error values according to the number of trees for each forest
#           --> the formatted dataframe used to plot
#           --> a list of prediction of the out-of-bag error for each forest (cross-val, only if cross_validation=T)


param.estimation.random.forest2<-function(dir="",Nindiv,intv_qtil=1e-4,nval=1000,n_tree=500,CI=c(0.025,0.975),cross_validation=T,paral=F,densityPlot=T){
  
  require(abcrf)
  # We need the functions from stefano to compute mpd + TD. 

  mpd_from_sfs<-function(sfs, folded=TRUE){
    ### check if it is folded or not
    if (isTRUE(folded)){n_ind<-2*length(sfs)}
    else {n_ind<-length(sfs)}
    ###
    mpd=0
    for (i in 1:length(sfs)){
      mpd=mpd+(((n_ind-i)*i)*sfs[i])
    }
    mpd=mpd/((n_ind*(n_ind-1))/2)
  
    return(mpd)
  }

  calcola_TD_folded<-function(folded_sfs){
    sample_size<-2*length(folded_sfs)
    theta_P<-mpd_from_sfs(folded_sfs)
    S<-sum(folded_sfs)

    ###### 
    calcola_a1<-function(sample_size){
        a<-0
        for (i in 1:(sample_size-1)){
            a<-a+(1/i)
        }
        return(a)
    }
    ##################

    ###### 
    calcola_a2<-function(sample_size){
        a<-0
        for (i in 1:(sample_size-1)){
            a<-a+(1/(i*i))
        }
    return(a)
    }
    ##################


    ###### 
    calcola_b1<-function(sample_size){

        a<-(sample_size+1)/(3*(sample_size-1))

        return(a)
    }
    ##################

    ###### 
    calcola_b2<-function(sample_size){

        a<-(2*((sample_size*sample_size)+sample_size+3))/(9*sample_size*(sample_size-1))

        return(a)
    }
    ##################

    ###### 
    calcola_c1<-function(sample_size){
        a<-calcola_b1(sample_size)-(1/calcola_a1(sample_size))
        return(a)
    }
    ##################


    ###### 
    calcola_c2<-function(sample_size){
        a<-calcola_b2(sample_size)-((sample_size+2)/(calcola_a1(sample_size)*sample_size))+(calcola_a2(sample_size)/(calcola_a1(sample_size)^2))

        return(a)
    }
    ##################

    calcola_e1<-function(sample_size){
        a<-calcola_c1(sample_size)/calcola_a1(sample_size)
        return(a)
    }

    calcola_e2<-function(sample_size){
        a<-calcola_c2(sample_size)/((calcola_a1(sample_size)^2)+calcola_a2(sample_size))
        return(a)
        }

    TD<-(theta_P-(S/calcola_a1(sample_size)))/(((calcola_e1(sample_size)*S)+(calcola_e2(sample_size)*S*(S-1)))^0.5)
    return(TD)
  }

  cat("==> Loading and computing sumstats...","\n")
  priors<-read.table(paste0(dir,"/priors.txt"),header=TRUE)
  sumstat<-read.table(paste0(dir,"/sumstat.txt"))
  if (ncol(sumstat)!=Nindiv){stop(paste0("sumstat file doesn't correspond to ",Nindiv," individuals folded-SFS"))}
  sumstat<-cbind(sumstat,apply(sumstat,1,mpd_from_sfs),apply(sumstat,1,sum),apply(sumstat,1,calcola_TD_folded))
  target<-read.table(paste0(dir,"/target.txt"))
  if (ncol(target)!=Nindiv){stop(paste0("target file doesn't correspond to a ",Nindiv," individuals folded-SFS"))}
  target<-cbind(target,apply(target,1,mpd_from_sfs),apply(target,1,sum),apply(target,1,calcola_TD_folded))
  cat("     Done.","\n")



  q<-seq(0, 1, by = intv_qtil)
  q<-as.vector(q)
  
  cross_val<-list()
  error<-list()
  data_list<-list()
  prediction_list<-list()
  stat_table<-c()
  if (cross_validation==F){
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
    pdf(paste0(dir,"/OOB_",colnames(priors)[i],".pdf"))
    plot(error[[i]])
    dev.off()
    cat("     Done.","\n","\n")
    
    if (cross_validation==T){
      cat("     >> computing cross-validation...","\n")
      cross_val[[i]]<-predictOOB(rf,training = data_list[[i]],paral=paral)
      cat("     Done.","\n","\n")
      
      k
    }
    if (densityPlot==T){
      cat("     >> printing abcrf density plot...","\n")
      pdf(paste0(dir,"/densityPlot_abcrf_",colnames(priors)[i],".pdf"))
      densityPlot(rf,obs = target,training = data_list[[i]],paral=paral)
      dev.off()
      cat("     Done.","\n","\n")
    }
    if (cross_validation==F){
      if (save_rf==T){
        rf_list[[i]]<-rf;names(rf_list)[i]<-colnames(priors)[i]
      } 
    }
    rm(rf)
  }
  if (cross_validation==T){names(cross_val)=colnames(priors)}
  
  names(prediction_list)=names(data_list)=names(error)=colnames(priors)
  
  df<-formatting.RF.est(prediction_list,dir = dir,byqtil=intv_qtil,nval=nval,CI=NULL)
  
  estimations<-estimate.plot.RF(df,priors,save = T,dir = dir,titles_param_plot = NULL)
  
  # For confidence intervalle
  CInames<-paste0(as.character(CI*100),"%")
  CI<-c(which(df$Quantiles>=CI[1])[1],which(df$Quantiles>=CI[2])[1])
  if(is.na(CI[2])){CI[2]<-nrow(df)}
  ic<-df[CI,]
  ic<-t(ic)
  colnames(ic)<-c("2.5%","97.5%")
  ic<-ic[-nrow(ic),]
  
  mediane<-c()
  mode<-c()
  for (i in 1:length(prediction_list)){
    mediane<-c(mediane,prediction_list[[i]]$med)
    temp<-density(df[,i]);mode_temp<-temp$x[which(temp$y==max(temp$y))]
    mode<-c(mode,mode_temp)
  }
  
  est<-data.frame(ic[,1],mode,mediane,ic[,2]);colnames(est)<-c(CInames[1],"Mode","Median",CInames[2])
  
  write.table(file = paste0(dir,"/estimates.txt"),est)
  
  #e<-err.regAbcrf(r.rf_fim_cuv_c, Nm_fim_cuv_c, paral = T)
  
  if (cross_validation==F) {
    fin<-list(data_list,rf_list,est,error,df);names(fin)<-c("ref_table","random_forests","estimated_values","error","reg_df")
  } else{fin<-list(data_list,est,error,cross_val,df);names(fin)<-c("ref_table","estimated_values","error","cross_validation","reg_df")}
  
  
  return(fin)
  
}