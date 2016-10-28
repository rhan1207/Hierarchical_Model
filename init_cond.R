init_cond <- function(mu_0, Lambda_0, eta_0, S_0, v_0, sigma_sqr_0){
  mu = mvrnorm(n = 1, mu_0, Lambda_0, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
  Sigma = riwish(eta_0, S_0)
  #alpha: shape parameter, beta: scalar rate parameter
  sigma_sqr = rinvgamma(1, v_0/2, v_0*sigma_sqr_0/2)     
  
  res = list(mu, Sigma,sigma_sqr)
  names(res) = c("mu", "Sigma", "sigma_sqr")
  return(res)
}