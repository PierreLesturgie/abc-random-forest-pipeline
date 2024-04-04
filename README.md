# ABC with Random Forest pipeline 
--- 

### R pipeline to run (and visualize) model selection and parameters estimation using abcRF package, with additional custom cross validation. 
####  Author: Pierre Lesturgie. Please contact (pierrelesturgie@outlook.fr) before using and sharing. 

This pipeline was first used and developped for the following (but updated in 04/2024): 

Lesturgie, P., Planes, S., & Mona, S. (2022). Coalescence times, life history traits and conservation concerns: An example from four coastal shark species from the Indo‐Pacific. *Molecular Ecology Resources*, 22(2), 554–566. https://doi.org/10.1111/1755-0998.13487

Lesturgie, P., Lainé, H., Suwalski, A., Chifflet-Belle, P., Maisano Delser, P., Clua, E., Jaquemet, S., Magalon, H., & Mona, S. (2022). Ecological and biogeographic features shaped the complex evolutionary history of an iconic apex predator (*Galeocerdo cuvier*). *BMC Ecology and Evolution*, 22(1), 147. https://doi.org/10.1186/s12862-022-02100-y

Lesturgie, P., Braun, C. D., Clua, E., Mourier, J., Thorrold, S. R., Vignaud, T., Planes, S., & Mona, S. (2023). Like a rolling stone: Colonization and migration dynamics of the gray reef shark (*Carcharhinus amblyrhynchos*). *Ecology and Evolution*, 13(1), 1–15. https://doi.org/10.1002/ece3.9746

---


## (1) Model Selection
### Requirements:
- functions from Stefano to compute mpd and TD
- A folder with sumstats and target

```
model.selection.random.forest(models=c("FIM","SST","NS"),nind,directory=c(""),analysis=c("all"),ntree=500,predictions=T,compute_oob=T,importance_variable=T,subset=NULL,LDA_plot=T)
```

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



## (2) Sumstat PCA
```
pcabc(models=c("FIM","SST","NS"),nind,directory=c(""),subset=NULL,pcs=NULL)
```
### Input : 
- models : names of models to compute (only important for loading sumstats and for output files)
- dir : directory for INPUT and output files. Must include **sumstat** files named *sumstat'model'.txt* (e.g: sumstat_FIM.txt) ***AND*** target summary statistics named **target.txt**
- nind : number of individuals
- subset : number of observation to subset each dataframe.
- pcs : number of pcs (if NULL will be asked) 
### Output : 
- write pdf of the pca
- returns a pca object


