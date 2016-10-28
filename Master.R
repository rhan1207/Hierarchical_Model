########################################################################
### Apply Gibbs Sampler to estimate parameters of Hierarchical Model ###
########################################################################
rm(list = ls())

library(MASS)
library(MCMCpack)
#root = "C:/Users/b1rxh09/Desktop/Hierarchical_Model"
root = "/Users/RJ/Desktop/Heirarchical_Model"
setwd(root)
## Read data and output ols coeffs and variance of target variable
source("process_data.R")
source("init_cond.R")


folder = "2015"
data_dir = paste(root, "Data", folder, sep="/")

setwd(data_dir)
all_files = list.files()

#create a matrix for all features' beta
#days = length(all_files)
days = 100
beta = matrix(0, 15, days)
tot_var = 0
for (t in (1:days)){ 
  file = all_files[t]
  print(file)
  temp = process_data(file,data_dir)  
  beta[,t] = temp$beta
  tot_var = tot_var + temp$var
}

# set hyper-parameters
ave_var = tot_var / days
cov_beta = cov(t(beta))
mean_beta = rowMeans(beta)
p = dim(beta)[1]
#initial conditions
theta = init_cond(mean_beta, cov_beta, p+2, cov_beta, 1, ave_var)

setwd(root)

