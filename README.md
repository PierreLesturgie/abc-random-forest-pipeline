# ABC with Random Forest pipeline 
--- 

####  author  = Pierre Lesturgie ****** ##
#### version = 0.1 ****** ##
####  email = pierrelesturgie@outlook.fr ****** ##
####  date = 11.02.2020 ****** ##

---

This is a R pipeline to run (and visualize) model selection and parameters estimation using abcRF package, with additional custom cross validation. 


#-------function 1 ==> performs the model selection

## (1) : Model Selection
### Requirements:
- functions from Stefano to compute mpd and TD
- A folder with sumstats and target

    missing
model.selection.random.forest<(models=c("FIM","SST","NS"),nind,directory=c(""),analysis=c("all"),ntree=500,predictions=T,compute_oob=T,importance_variable=T,subset=NULL,LDA_plot=T)


### Input : 
- models : names of models to compute (only important for loading sumstats and for output files)
- dir : directory for INPUT and output files. Must include **sumstat** files named *sumstat'model'.txt* (e.g: sumstat_FIM.txt) ***AND*** target summary statistics named **target.txt**
- nind : number of individuals
- analysis : "LDA", "basic" or "all", to compute only the model selection with a linear discriminant analysis, without, or both.
- ntree : number of trees used to grow the forest
- predictions : does the function computes model choice (or only grows forest) ?
- compute_oob : does the function computes the evolution of out-of-bag error with the number of trees ?
- importance_variable : does the function tests the importance of the variables in the decision ?
- subset : number of observation to subset each dataframe.
- LDA_plot : does the function computes the LDA plot (very similar to a pca)
### Output : 
- returns a list with the confusion matrix and the output of analyses aseked to the function. 



## (2) : Sumstat PCA

pcabc(models=c("FIM","SST","NS"),nind,directory=c(""),subset=NULL,pcs=NULL)

### Input : 
- models : names of models to compute (only important for loading sumstats and for output files)
- dir : directory for INPUT and output files. Must include **sumstat** files named *sumstat'model'.txt* (e.g: sumstat_FIM.txt) ***AND*** target summary statistics named **target.txt**
- nind : number of individuals
- subset : number of observation to subset each dataframe.
- pcs : number of pcs (if NULL will be asked) 
### Output : 
- write pdf of the pca
- returns a pca object


