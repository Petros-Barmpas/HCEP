# bicUMAP
Biscecting Hierarchical Clustering 

## Introduction

The ATHLOS cohort is composed of several harmonized datasets of international cohorts related to health and aging. 
The healthy aging index has been constructed based on a selection of particular variables from 16 individual studies. 
In this paper, we consider a selection of additional variables found in ATHLOS and investigate their utilization for 
predicting the healthy aging index. For this purpose motivated by the volume and diversity of the dataset, we focus 
our attention upon the clustering for prediction scheme, where unsupervised learning is utilized to enhance prediction 
power, showing the predictive utility of exploiting structure in the data by clustering. We show that imposed computation 
bottlenecks can be surpassed when using appropriate hierarchical clustering, within a clustering for ensemble classification 
scheme, while retaining prediction benefits. We propose a complete methodology that is evaluated against baseline methods 
and the original concept. The results are very encouraging suggesting further developments in this direction along with 
applications in tasks with similar characteristics. A straightforward open source implementation for the R project is provided.

![] https://github.com/Petros-Barmpas/bicUMAP/blob/master/output.png
*Clustering example procedure*

## Reference
If you use this code please reference the corresponding recently published article. 

Barmpas, Petros & Tasoulis, Sotiris & Vrahatis, Aristidis & Prina, Matthew & Ayuso-Mateos, José & Bickenbach, Jerome & Bayés, Ivet & Bobak, Martin & Caballero, Francisco Félix & Chatterji,
 Somnath & Egea-Cortés, Laia & Garcia-Esquinas, Esther & Leonardi, Matilde & Koskinen, Seppo & Koupil, Ilona & Pajak, Andrzej & Prince, Martin & Sanderson, Warren & Scherbov, Sergei & Panagiotakos, 
Demosthenes. (2021). A Hybrid Machine Learning Framework for Enhancing the Prediction Power in Large Scale Population Studies: The ATHLOS Project. 10.1101/2021.01.23.21250355. 

[download bibtex](https://github.com/Petros-Barmpas/bicUMAP/blob/master/bibtex.txt)

## License
This project is licensed under the BSD-3-Clause License - see the LICENSE file for details.

## Example script
```r
# Example Script


# CRIM - per capita crime rate by town
# ZN - proportion of residential land zoned for lots over 25,000 sq.ft.
# INDUS - proportion of non-retail business acres per town.
# CHAS - Charles River dummy variable (1 if tract bounds river; 0 otherwise)
# NOX - nitric oxides concentration (parts per 10 million)
# RM - average number of rooms per dwelling
# AGE - proportion of owner-occupied units built prior to 1940
# DIS - weighted distances to five Boston employment centres
# RAD - index of accessibility to radial highways
# TAX - full-value property-tax rate per $10,000
# PTRATIO - pupil-teacher ratio by town
# B - 1000(Bk - 0.63)^2 where Bk is the proportion of blacks by town
# LSTAT - % lower status of the population
# MEDV - Median value of owner-occupied homes in $1000's


install.packages("mlbench")

# source the file including all required functions for dePDDP


source("dePDDP.R")
source("ensemble_functions.R")

# load required libraries

library(mlbench)
library(knitr)
library(randomForest)
library(doMC)
library(BBmisc)
library(irlba)

# load the UCI’s Boston House Dataset

data(BostonHousing)

#became familliar with the dataset
dim(BostonHousing)
head(BostonHousing)


# look at the pairs to better understand the relationship between variables
pairs(BostonHousing, col='red')


# look into the correlations between the numeric variables and show them
housing_data_nc = BostonHousing[, -4]

corrmatrix = cor(housing_data_nc)

kable(t(corrmatrix))


## setting the seed to make the example reproducible
set.seed(1234)


# Visuilize the dataset in two dimensions to estimate the number of clusters in this case
embedumap <- uwot:: umap(BostonHousing, n_neighbors = 50, min_dist = 0.3, n_components = 2)
embedumap <- as.data.frame(embedumap)

cl_n <- 3

# Visuilize the dataset in two dimensions to estimate the number of clusters in this case
embedumap <- uwot:: umap(BostonHousing, n_neighbors = 50, min_dist = 0.3, n_components = 5)
embedumap <- as.data.frame(embedumap)

# res<-dePDDP(BostonHousing, k=cl_n,SP="DEN")

smp_size <- floor(0.25 * nrow(BostonHousing))


#shuffle the data
workData <- shuffle_Data(BostonHousing, seed = 1234)

workClass <- workData$shuffled_class

workData <- workData$shuffled_data

#transofming the factor variable as boolean integer
workData$chas <- as.integer(as.character(workData$chas))

# cluster the dataset ===
# ClustID_Export <- onetime_clustering_export_instances(workData,kk = cl_n, tsamp_size = smp_size)
ClustID_Export <- onetime_clustering_export_instances(embedumap,kk = cl_n, tsamp_size = smp_size)


# Linear Regression ensemble
LregResult <- linRegEnsemble(workData,idExp = ClustID_Export,respVar = "medv")
print(paste("The resulted RMSE score is:" , LregResult$rmses_mean_asc[length(LregResult$rmses_mean_asc)]))
print(paste("The resulted RSQ score is:" , LregResult$rsq_mean_asc[length(LregResult$rsq_mean_asc)]))

#Random Forest Regression
RFregResult <- RFregEnsemble(workData,idExp = ClustID_Export,respVar = "medv")

print(paste("The resulted RMSE score is:" , RFregResult$rmses_mean_asc[length(RFregResult$rmses_mean_asc)]))
print(paste("The resulted RSQ score is:" , RFregResult$rsq_mean_asc[length(RFregResult$rsq_mean_asc)]))

#


```

