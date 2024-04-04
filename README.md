# ABC with Random Forest pipeline 
--- 

### R pipeline to run (and visualize) model selection and parameters estimation using abcRF package, with additional custom cross validation. 
####  Author: Pierre Lesturgie. Please contact (pierrelesturgie@outlook.fr) before using and sharing. 

This pipeline was first used and developped for the following (but updated in 04/2024): 

Lesturgie, P., Planes, S., & Mona, S. (2022). Coalescence times, life history traits and conservation concerns: An example from four coastal shark species from the Indo‐Pacific. *Molecular Ecology Resources*, 22(2), 554–566. https://doi.org/10.1111/1755-0998.13487

Lesturgie, P., Lainé, H., Suwalski, A., Chifflet-Belle, P., Maisano Delser, P., Clua, E., Jaquemet, S., Magalon, H., & Mona, S. (2022). Ecological and biogeographic features shaped the complex evolutionary history of an iconic apex predator (*Galeocerdo cuvier*). *BMC Ecology and Evolution*, 22(1), 147. https://doi.org/10.1186/s12862-022-02100-y

Lesturgie, P., Braun, C. D., Clua, E., Mourier, J., Thorrold, S. R., Vignaud, T., Planes, S., & Mona, S. (2023). Like a rolling stone: Colonization and migration dynamics of the gray reef shark (*Carcharhinus amblyrhynchos*). *Ecology and Evolution*, 13(1), 1–15. https://doi.org/10.1002/ece3.9746

---

The pipeline performs model selection and parameter estimation using an Approximate Bayesian Computation with random forests approach. 
It consists in two steps: 

1. Model selection: estimating which model best fits the observed data (using a classification algorithm)
2. Parameter estimation: parameter estimates (mean, median, mode, CI) under a specific model (using a regression algorithm). 

A script to run an example and according data is present in the folder "test"

## (1) Model Selection

Performs  model selection from a target + a named list of summary statistics simulated from different scenarios. 

```

ms<-model.selection.random.forest(target = target,list_sumstat = sumstat_list,
                                           directory = "./",
                                           analysis = "all",ntree = 500,predictions = T,compute_oob = F,
                                           importance_variable = F,subset=NULL,LDA_plot=T,paral = F)
```

### List of arguments: 
- target : observed summary statistics (e.g., SFS, genetic diversity, Tajima's D,...)
- list_sumstat : Named list of summary statistics obtained from simulated models 
- dir : directory output files.
- analysis : "LDA", "basic" or "all", to compute the ms adding linear discriminant analysis axes as sumstat, without, or both. 
- ntree : number of trees used to grow the forest
- predictions : does the function computes model choice (or only grows forest) ?
- compute_oob : does the function computes the evolution of out-of-bag error with the number of trees ?
- importance_variable : does the function tests the importance of the variables in the decision ?
- subset : number of observation to subset each dataframe.
- LDA_plot : does the function computes the LDA plot (very similar to a pca)
- Paral : to parallelize using ncores - 1 

The function returns a list with the confusion matrix and the output of analyses aseked to the function. 

### Optionally, one can run a PCA with the target and the sumstat, to check consistency between them 
However, it is also recommended to do the LDA plot in the model.selection.random.forest() in place of the PCA 

```
pca_sumstats = pcabc(target = target, list_sumstat = sumstat_list,directory = "./",pcs = NULL)
```

#### List of arguments: 
- target : observed summary statistics (e.g., SFS, genetic diversity, Tajima's D,...)
- list_sumstat : Named list of summary statistics obtained from simulated models 
- dir : directory output files.
- subset : number of observation to subset each dataframe.
- pcs : number of pcs (if NULL will be asked) 

Returns a pca object and writes pdf of the pca

## (2) Parameter estimation 

```
est = param.estimation.random.forest(sumstat = sumstat,priors = priors,target = target, 
                               dir = "./",n_tree = 500,
                               cross_validation = T,npods=100,paral=F,write_cv=F)
```

Performs parameter estimation, and optionally cross-validation

### List of arguments: 
- sumstat : data frame w/ sumstat data (SFS + other wanted sumstats)
- priors : priors corresponding to the sumstat data (i.e., the simulated parameters, such as effective size, migration rate, etc.)
- target : targetted sumstat (need to have exactly the same stats computed than sumstat)
- intv_qtil : steps between quantiles (for the plots)
- nval : number of points used to plot the values from the quantiles
- n_tree : number of trees used to grow the forest
- CI : bounds of the confidence interval
- dir : directory for output files
- cross_validation : TRUE for computation of cross-validation
- npods = number of Pseudo-Observed-Datasets sampled for cross-validation
- predict_OOB_CV : computation of abcRF package implemented cross validation <<< very time consuming >>>
- densityPlot : TRUE to compute and save the abcRF implemented density plot
- save_rf : TRUE to save the random forests (only if cross_validation=F) <<< very heavy objects >>>

The function returns a list with *reference tables*; *estimated values*; *out-of-bag* error values; formated *dataframe* for plots; *cross-validation* results; and the *random_forest* . Additionally, it writes results and prints graphs in the **dir**.



