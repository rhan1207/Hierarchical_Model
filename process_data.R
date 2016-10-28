#Inputs: daily stock return with Beta, momemtum, size, volotility and value factors
#Outputs: linear regression coefficients and variance for stock return
process_data <- function(file,data_dir){
  
  setwd(data_dir)
  #df <- read.table("apt.2015-06-30.txt", header = TRUE, sep = "|")
  df <- read.table(file, header = TRUE, sep = "|")

  # delete row without industry
  del = which(df$INDUSTRY == "Unknown")
  if (length(del) > 0){df = df[-del,]} 
  # create dummy variables for all industries
  for(indstr in unique(df$INDUSTRY)){
    df[paste("INDUSTRY", indstr, sep="_")] <- ifelse(df$INDUSTRY == indstr, 1, 0)  
  }
  
  df <- df[,-c(1:2, 8:11)]
  
  # create compute ols coefficients
  beta <- lm(df$RETNEXT ~ 0 +., data = df)
  beta <- beta$coefficients
  # calculate variance of y
  var_ret = var(df$RETNEXT)

  out = list(beta, var_ret)
  names(out) = c("beta", "var")
  
  return(out)
  
  
}