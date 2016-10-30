Gibbs_Sampler <- function(theta, data, iter = 10000){
  
  m = theta$days
  p = length(theta$mu)
  
  Beta_mat = matrix(0, nrow = p, ncol = m)
  
  Mu_mat = matrix(0, nrow = iter, ncol = p)
  
  # only the diagonal elements of Sigma
  Sigma_mat = matrix(0, nrow = iter, ncol = p)
  
  # variance for y
  sigma_vec = matrix(0, nrow = iter, ncol = 1)
  
  
  # initial values for Gibbs Sampling
  Sigma = theta$Sigma
  mu = matrix(theta$mu, p, 1)
  sigma_sqr = theta$sigma_sqr
  mu_0 = theta$mu_0
  S_0 = theta$S_0
  eta_0 = theta$eta_0
  v_0 = theta$v_0
  sigma_sqr_0 = theta$sigma_sqr_0
  Lambda_inv = solve(theta$Lambda)
  
  
  for (t in (1:iter)){
    print(t)
    #during the first pass, cache (X'X)^(-1)X'Y
    X_Y_Mat = matrix(0, p, m)
    
    # cache these values for sampling all beta
    Sigma_mu = solve(Sigma)%*%matrix(mu, p, 1)
    sigma_sqr_inv = 1/sigma_sqr * diag(p)
    Sigma_inv = solve(Sigma)
    
    # compute total obs and SSR for sampling sigma^2 
    tot_obs = 0
    SSR = 0
    ## Sample beta ##
    for (i in 1:m){
      if (t == 1){
        Y = data[[i]][,6, drop = FALSE]
        X = data[[i]][,-6]
        X_Y = solve(t(X)%*%X)%*%t(X)%*%Y
        X_Y_Mat[,i] = X_Y
        tot_obs = tot_obs + dim(Y)[1]
      }
      
      #cache Sigma', mu*Sigma^(-1) and sigma^{-2}*I_p
      Sigma_1 = solve(Sigma_inv + sigma_sqr_inv)
      mu_1 = Sigma_1 %*% (Sigma_inv%*%mu + sigma_sqr_inv %*% X_Y_Mat[,i, drop = FALSE])
      beta = mvrnorm(n = 1, mu_1, Sigma_1, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
  
      Beta_mat[,i] = beta
      SSR = SSR + sum((Y - X %*% Beta_mat[,i])*(Y - X %*% Beta_mat[,i]))
    }
    
    
    ## Sample mu ##
    Sigma_p = solve(Lambda_inv + m * Sigma_inv)
    mu_p = Sigma_p%*%(Sigma_inv%*%rowSums(Beta_mat) + Lambda_inv %*% matrix(mu_0, p, 1))
    mu = mvrnorm(n = 1, mu_p, Sigma_p, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
    
    mu = matrix(mu, p, 1)
    Mu_mat[t,] = mu
    
    
    
    ## Sample Sigma ##
    S = S_0
    for (j in (1 : m)){
      S = S + (mu - Beta_mat[,j,drop = FALSE])%*%t(mu - Beta_mat[,j,drop = FALSE])
    }
    Sigma = riwish(eta_0 + m, S)
    Sigma_mat[t, ] = diag(Sigma)
    
    
    # Sample sigma^2 ##
    sigma_sqr = rinvgamma(1, v_0/2 + tot_obs/2, (v_0*sigma_sqr_0 + SSR)/2)
    sigma_vec[t, 1] = sigma_sqr
  }
  
  res = list(Mu_mat, Sigma_mat, sigma_vec)
  return(res)
}