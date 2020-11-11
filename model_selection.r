################################################################
#------- FUNCTIONS MODEL SELECTION W/ RANDOM FORESTS ABC ------#
################################################################


#-------function 1 ==> performs the model selection

# >>>>>> function 1 : performs the model selection

#       <<< Needs functions from Stefano to compute mpd and TD >>>
#       <<< Needs a folder with sumstats and target >>>

# ==> Input : 
#       *** models : names of models to compute (only important for loading sumstats and for output files)
#       *** dir : directory for INPUT and output files
#           --> sumstat files named "sumstat'model'.txt" (ex sumstat_FIM.txt) << only the folded-SFS >>
#           --> target file named "target.txt" << only the folded-SFS >>
#       *** nind : number of individuals
#       *** analysis : "LDA", "basic" or "all", to compute only the ms with a linear discriminant analysis, without, or both. 
#       *** ntree : number of trees used to grow the forest
#       *** predictions : does the function computes model choice (or only grows forest) ?
#       *** compute_oob : does the function computes the evolution of out-of-bag error with the number of trees ? 
#       *** importance_variable : does the function tests the importance of the variables in the decision ? 
#       *** subset : number of observation to subset each dataframe. 
#       *** LDA_plot : does the function computes the LDA plot (very similar to a pca)
# ==> Output : 
#       *** returns a list with the confusion matrix and the output of analyses aseked to the function. 


model.selection.random.forest<-function(models=c("FIM","SST","NS"),nind,directory=c(""),analysis=c("all"),ntree=500,predictions=T,compute_oob=T,importance_variable=T,subset=NULL,LDA_plot=T){
  require("abcrf")
  
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
  
  target<-read.table(paste0(directory,"/target.txt"))
  if (length(target)!=nind){stop(paste0("The target doesn't correspond to a ",nind," individuals folded SFS"))}
  target<-as.double(cbind(target, mpd_from_sfs(target),sum(target),calcola_TD_folded(target)))
  target<-as.data.frame(t(target))
  
  sumstat<-c()
  nsims<-c()
  for (i in 1:length(models)){
    temp<-read.table(paste0(directory,"/sumstat_",models[i],".txt"))[,c(1:nind)]
    if(is.null(subset)==F){
      m<-c(sample(x = c(1:nrow(temp)),size = subset,replace = F))
      temp<-temp[m,]
    }
    sumstat<-rbind(sumstat,temp)
    nsims<-rbind(nsims,nrow(temp))
  }
  
  
  cat("Computing summary statistics on the simulated dataset...")
  MPD<-apply(sumstat,1,mpd_from_sfs)
  S<-apply(sumstat,1,sum)
  TD<-apply(sumstat,1,calcola_TD_folded)
  sumstat<-cbind(sumstat,MPD,S,TD)
  cat(" Done.","\n")
  
  if (length(which(is.na(sumstat$TD)))>=1){
    sumstat[which(is.na(sumstat$TD)),"TD"]<-0
  }
  
  y<-c()
  for (i in 1:length(models)){
    y<-rbind(y,as.matrix(rep(models[i],nsims[i])))
  }
  y<-as.vector(y)
  
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
      pLDA<-predict(rfLDA, target, training = data_rf, ncores = 10, ncores.predict = 10,ntree = ntree)
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
      p<-predict(rf, target, training = data_rf, ncores = 10, ncores.predict = 10,ntree = ntree)
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
      pLDA<-predict(rfLDA, target, training = data_rf, ncores = 10, ncores.predict = 10,ntree = ntree)
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
      p<-predict(rf, target, training = data_rf, ncores = 10, ncores.predict = 10,ntree = ntree)
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
    y2<-as.vector(c(y,"TARGET"))
    sumstat2<-rbind(sumstat,target)
    data_rf2<-data.frame(y2,sumstat2)
    rf2<-abcrf(y2~., data = data_rf2, 
               lda=T, ntree=ntree,
               paral=FALSE)
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


# >>>>>> function 1 : performs a pca for the abc

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

pcabc<-function(models=c("FIM","SST","NS"),nind,directory=c(""),subset=NULL,pcs=NULL){
  require(ade4)
  target<-read.table(paste0(directory,"/target.txt"))
  if (length(target)!=nind){stop(paste0("The target doesn't correspond to a ",nind," individuals folded SFS"))}
  target<-as.double(cbind(target, mpd_from_sfs(target),sum(target),calcola_TD_folded(target)))
  target<-as.data.frame(t(target))
  
  sumstat<-c()
  nsims<-c()
  for (i in 1:length(models)){
    temp<-read.table(paste0(directory,"/sumstat_",models[i],".txt"))[,c(1:nind)]
    if(is.null(subset)==F){
      m<-c(sample(x = c(1:nrow(temp)),size = subset,replace = F))
      temp<-temp[m,]
    }
    sumstat<-rbind(sumstat,temp)
    nsims<-rbind(nsims,nrow(temp))
  }
  
  
  cat("Computing summary statistics on the simulated dataset...")
  MPD<-apply(sumstat,1,mpd_from_sfs)
  S<-apply(sumstat,1,sum)
  TD<-apply(sumstat,1,calcola_TD_folded)
  sumstat<-cbind(sumstat,MPD,S,TD)
  cat(" Done.","\n")
  
  if (length(which(is.na(sumstat$TD)))>=1){
    sumstat[which(is.na(sumstat$TD)),"TD"]<-0
  }
  
  y<-c()
  for (i in 1:length(models)){
    y<-rbind(y,as.matrix(rep(models[i],nsims[i])))
  }
  y<-as.vector(y)
  
  data_rf<-data.frame(y,sumstat)
  colnames(target)<-colnames(data_rf)[-1]
  
  y<-as.vector(c(y,"TARGET"))
  sumstat<-rbind(sumstat,target)
  data_rf<-data.frame(y,sumstat)
  
  if (is.null(pcs)==F){
    pca<-dudi.pca(df = data_rf[, -1], scannf = FALSE, nf = pcs)
  }
  pca<-dudi.pca(data_rf[,-1])
  
  col<-c()
  for (i in 1:length(models)){
    col=c(col,rep(i,nsims[i]))
  }
  
  pcs<-pca$nf
  d<-(pca$eig/sum(pca$eig))*100
  
  pdf(paste0(directory,"pca.pdf"),onefile = T)
  for (j in 1:(pcs-1)){
    plot(pca$tab[-1,c(1,j+1)],col=col,pch = 19,cex=0.2,xlab=paste0("PC1 (",round(d[1]),"%)"),
         ylab=paste0("PC",j+1," (",round(d[j+1]),"%)")) + points(pca$tab[1,c(1,j+1)],
                                                                 col=i+1,pch = 18,cex=1)
    legend("topleft",legend = c(models,"Target"),col=c(1:(i+1)),pch = 19)
    
  }
  
  dev.off()
  return(pca)
}
