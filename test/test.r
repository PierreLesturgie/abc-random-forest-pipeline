library(abcrf)
source("../abc_pipeline.r")


target = read.table("target.txt")

sumstat_list=list(
  read.table("sumstat_NS.txt"),
  read.table("sumstat_FIM.txt"),
  read.table("sumstat_SST.txt")
)
names(sumstat_list) = c("NS","FIM","SST")

### can do the PCA, but otherwise recommended to do the LDA plot in the model.selection.random.forest()
pca_sumstats = pcabc(target = target, list_sumstat = sumstat_list,directory = "./",pcs = 2)

ms<-model.selection.random.forest(target = target,list_sumstat = sumstat_list,
                                           directory = "./",
                                           analysis = "all",ntree = 500,predictions = T,compute_oob = F,
                                           importance_variable = F,subset=NULL,LDA_plot=T,paral = F)


### SST wins 
priors = read.table("priors_SST.txt")
sumstat = sumstat_list$SST

est = param.estimation.random.forest(sumstat = sumstat,priors = priors,target = target, 
                               dir = "./",n_tree = 500,
                               cross_validation = T,npods=100,paral=F,write_cv=F)

