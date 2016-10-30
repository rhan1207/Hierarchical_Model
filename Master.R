########################################################################
### Apply Gibbs Sampler to estimate parameters of Hierarchical Model ###
########################################################################
rm(list = ls())

library(MASS)
library(MCMCpack)
root = "/Users/RJ/GitHub/Hierarchical_Model"
setwd(root)
## Read data and output ols coeffs and variance of target variable
source("process_data.R")
source("init_cond.R")
source("Gibbs_Sampler.R")

year = c(2006:2014)
for (y in year){
folder = toString(y)
data_dir = paste(root, "Data", folder, sep="/")

setwd(data_dir)
all_files = list.files()

#create a matrix for all features' beta
days = length(all_files)
beta = matrix(0, 15, days)
data = list()
tot_var = 0
for (t in (1:days)){ 
  file = all_files[t]
  print(file)
  temp = process_data(file,data_dir)  
  beta[,t] = temp$beta
  tot_var = tot_var + temp$var
  data[[t]] = temp$data
}


setwd(root)
# set hyper-parameters
ave_var = tot_var / days
cov_beta = cov(t(beta))
mean_beta = rowMeans(beta)
p = dim(beta)[1]

#initial conditions
theta = init_cond(mean_beta, cov_beta, p+2, cov_beta, 1, ave_var)
theta$days = days
theta$mu_0 = mean_beta
theta$Lambda = cov_beta
theta$eta_0 = p+2
theta$S_0 = cov_beta
theta$v_0 = 1
theta$sigma_sqr_0 = ave_var

#run Gibbs Sampler
results = Gibbs_Sampler(theta, data)

quant_95 <- function(x){return(quantile(x, c(0.025, .975)))}

mu_q = apply(results[[1]], 2, quant_95)
cov_q = apply(results[[2]], 2, quant_95)
sig_q = apply(results[[3]], 2, quant_95)


write.table(mu_q, paste(folder, "mu.csv", sep="_"))
write.table(results[[1]], paste(folder, "mu.csv", sep="_"), row.names = FALSE, col.names = FALSE, append = TRUE)

write.table(cov_q, paste(folder, "cov.csv", sep="_"))
write.table(results[[2]], paste(folder, "cov.csv", sep="_"), row.names = FALSE, col.names = FALSE, append = TRUE)

write.table(sig_q, paste(folder, "sigma.csv", sep="_"))
write.table(results[[3]], paste(folder, "sigma.csv", sep="_"), row.names = FALSE, col.names = FALSE, append = TRUE)
}