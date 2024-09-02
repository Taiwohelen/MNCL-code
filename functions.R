#setwd("")

############
sparsity_generator <- function(p, prop, a, b, c, d){
  upper = (p^2 - p)/2
  density = ceiling(prop * upper)
  x = rep(0, density)
  y <- runif(density, 0, b-a+d-c)
  for (i in 1:density){
    if( y[i] < (b-a) ){
      x[i] <- a + y[i]
    }else{
      x[i] <- c + y[i] - (b-a)
    }
  }
  cell = c(x, rep(0, upper-density))
  sampledcell = sample(cell)
  mat = matrix(0, nrow = p, ncol = p)
  c = 1
  for(i in 1:(p-1)){
    for(j in (i+1):p){
      mat[i,j] = mat[j,i] = sampledcell[c]
      c = c + 1
    }
  }
  diag(mat) = 1
  return(mat)
}

# X : n X p covariate matrix      ##this is the matrix of n p-variate observations

############

BSSC_pMoM_Mode_max <- function(X, qn, r, s, niter, nburn){
  
  n = nrow(X)  #no of responses
  p = ncol(X)   #no of variate on each response 
  S = 1/n*t(X)%*%X  #matrix of error of X, where X has mean 0 
  
  res = list()  ##list to hold the output of the iteration
  
  Omega.mat = array(0, dim = c(p, p, niter + nburn))
  
  
  # set initial value for Omega.mat using estimates from glasso
  
  library(glasso)
  Vs = var(X)
  a = glasso(Vs, rho=.055)
  aa = glasso(Vs,rho=.1, w.init=a$w, wi.init=a$wi)  ##  a$w is estimated covariance matrix
  glassoaa = aa$wi                      ## aa$wi is estimated inverse covariance matrix
  Omega.mat[,,1] = glassoaa       ## first sample of the precision matrix for initial value
  
  target = function(wjk, b, a, n){
    return(dnorm(x = wjk, mean = -b/a, sd = sqrt(1/(n*a)))*wjk^2)
  }
  endtime <- rep(0, niter+nburn)
  # MCMC sampling
  for(i in 2:(niter+nburn)){ #for each MCMC sample, recall, we already have the first sample as glasso estimate
    #cat(i,"th iteration is completed from i. . . . . .\n"
    starttime = Sys.time()
    
    Omega = Omega.mat[,,i-1]
    L = matrix(0,p,p)    ### initial value for L -- this is the latent variable describing that w_jk...
    ##  ...is active or inactive. Hence a matrix of 0 and 1.
    ##   It follows a Bernoulli distribution with prob p_jk
    for(j in 1:(p-1)){
      #cat(j,"th iteration is completed from j. . . . . .\n")
      for(k in (j+1):p){
        #cat(k,"th iteration is completed from k . . . . . .\n") 
        #for off diagonal elements of the matrix
        lambda = rgamma(1,shape = r + 1.5, rate = (0.5*Omega[j,k]*Omega[j,k]) + s)  #posterior of off diagonal parameter lambda
        
        a = S[j,j] + S[k,k] + lambda/n   ##fixed constant for each lambda_jk
        b = t(Omega[j,])%*%S[k,] + t(Omega[,k])%*%S[j,] - Omega[j,k]*(S[j,j] + S[k,k])   ##fixed constant for each w_jk
        
        ## using log/threshold
        logc = log(qn) - log(1 - qn) + (1.5*log(lambda/(n*a))) + (log(1 + (n*b^2/a))) + (n*b^2/(2*a))
        
        # sample L_jk. This is the posterior for latent variable of active/inactive omega_jk
        if(logc[1,1] < 10){ L[j,k] = L[k,j] = rbinom(n = 1, size = 1, prob = 1 - (1/(1+exp(logc))))
        } else { L[j,k] = L[k,j] = rbinom(n = 1, size = 1, prob = 1 - 1/(1+exp(10)))
        }
        
        
        # updating the ith sample of Omega_jk based on the sample of the latent variable
        
        if(L[j,k] == 0) { Omega[j,k] = Omega[k,j] = 0
        }else {   ###choose the w_jk that maximises the target density
          w1 = (sqrt((n*b)^2+(8*n*a)) - (n*b)) / (2*n*a)
          w2 = (- sqrt((n*b)^2+(8*n*a)) - (n*b)) / (2*n*a)
          d1 = target(w1, b, a, n)
          d2 = target(w2, b, a, n)
          if(d1>d2){Omega[j,k] = Omega[k,j] = w1
          } else {Omega[j,k] = Omega[k,j] = w2}
          
        }
      }
      
      
      #for diagonal elements of the matrix
      bj = t(Omega[j,])%*%S[j,] - Omega[j,j]*S[j,j]
      gamma <- rgamma(1,shape = r + 1, rate = Omega[j,j] + s)  #posterior of diagonal parameter gamma
      
      Omega[j,j] = (sqrt((gamma + n*bj)^2 + 4*(n^2)*S[j,j]) - (gamma + n*bj))/(2*n*S[j,j])
    }
    
    
    bp = t(Omega[p,])%*%S[p,] - Omega[p,p]*S[p,p]
    gamma <- rgamma(1,shape = r + 1, rate = Omega[p,p] + s)  #posterior of the pth diagonal parameter gamma
    
    Omega[p,p] = (sqrt((gamma + n*bp)^2 + 4*(n^2)*S[p,p]) - (gamma + n*bp))/(2*n*S[p,p])
    Omega.mat[,,i] = Omega  ##the full matrix of the one generated sample of the Omega
    
    endtime[i] = Sys.time() - starttime
    
    #if(i %% 1000 == 0)
    #cat(i,"th iteration is completed. . . . . .\n")
  }
  
  # end of (i in 2:(niter+nburn)) for loop
  
  
  res = list(Omega.mat = Omega.mat, Time = endtime)
  
  return(res)
  
}

###########

BSSC_normal <- function(X, qn, r, s, niter, nburn){
  
  n = nrow(X)  #no of responses
  p = ncol(X)   #no of variate on each response 
  S = 1/n*t(X)%*%X  #matrix of error of X, where X has mean 0 
  
  res = list()  ##list to hold the output of the iteration
  
  Omega.mat = array(0, dim = c(p, p, niter + nburn))
  
  
  # set initial value for Omega.mat using estimates from glasso
  
  library(glasso)
  Vs = var(X)
  a = glasso(Vs, rho=.055)
  aa = glasso(Vs,rho=.1, w.init=a$w, wi.init=a$wi)  ##  a$w is estimated covariance matrix
  glassoaa = aa$wi                      ## aa$wi is estimated inverse covariance matrix
  Omega.mat[,,1] = glassoaa       ## first sample of the precision matrix for initial value
  #Omega.mat[,,1] = diag(p)
  
  endtime <- rep(0, niter+nburn)
  # MCMC sampling
  for(i in 2:(niter+nburn)){ #for each MCMC sample, recall, we already have the first sample as glasso estimate
    
    starttime = Sys.time()
    
    Omega = Omega.mat[,,i-1]
    L = matrix(0,p,p)    ### initial value for L -- this is the latent variable describing that w_jk...
    ##  ...is active or inactive. Hence a matrix of 0 and 1.
    ##   It follows a Bernoulli distribution with prob p_jk
    for(j in 1:(p-1)){
      
      for(k in (j+1):p){
        
        #for off diagonal elements of the matrix
        lambda = rgamma(1,shape = r + 0.5, rate = (0.5*Omega[j,k]*Omega[j,k]) + s)  #posterior of off diagonal parameter lambda
        
        a = S[j,j] + S[k,k] + lambda/n   ##fixed constant for each lambda_jk
        b = t(Omega[j,])%*%S[k,] + t(Omega[,k])%*%S[j,] - Omega[j,k]*(S[j,j] + S[k,k])   ##fixed constant for each w_jk
        
        ## using log/threshold
        logc = log(qn) - log(1 - qn) + (0.5*log(lambda/(n*a))) + (n*b^2/(2*a))
        
        # sample L_jk. This is the posterior for latent variable of active/inactive omega_jk
        if(logc[1,1] < 10){ L[j,k] = L[k,j] = rbinom(n = 1, size = 1, prob = 1 - (1/(1+exp(logc))))
        } else { L[j,k] = L[k,j] = rbinom(n = 1, size = 1, prob = 1 - 1/(1+exp(10)))
        }
        
        # updating the ith sample of Omega_jk based on the sample of the latent variable
        if(L[j,k] == 0) { Omega[j,k] = Omega[k,j] = 0
        }else {
          Omega[j,k] = Omega[k,j] = rnorm(n = 1, mean = -b/a, sd = sqrt(1/(n*a)))
        }
      }
      
      #for diagonal elements of the matrix
      bj = t(Omega[j,])%*%S[j,] - Omega[j,j]*S[j,j]
      gamma <- rgamma(1,shape = r + 1, rate = Omega[j,j] + s)  #posterior of diagonal parameter gamma
      
      Omega[j,j] = (sqrt((gamma + n*bj)^2 + 4*(n^2)*S[j,j]) - (gamma + n*bj))/(2*n*S[j,j])
    }
    
    
    bp = t(Omega[p,])%*%S[p,] - Omega[p,p]*S[p,p]
    gamma <- rgamma(1,shape = r + 1, rate = Omega[p,p] + s)  #posterior of the pth diagonal parameter gamma
    
    Omega[p,p] = (sqrt((gamma + n*bp)^2 + 4*(n^2)*S[p,p]) - (gamma + n*bp))/(2*n*S[p,p])
    Omega.mat[,,i] = Omega  ##the full matrix of the one generated sample of the Omega
    
    endtime[i] = Sys.time() - starttime
    
    #if(i %% 1000 == 0)
    # cat(i,"th iteration is completed. . . . . .\n")
  }
  
  # end of (i in 2:(niter+nburn)) for loop
  
  
  res = list(Omega.mat = Omega.mat, Time = endtime)
  
  return(res)
  
}

######### ESTIMATE OF MAGNITUDE / RE ############

BSSC_estimate = function(p, D, Omega_iter_result, Omega_iter_sparsity, final_Omega_sparsity){
  #Omega_iter_result =BSSC_nonlocal_Mode
  #Omega_iter_sparsity = BSSC_nonlocal_Omega_1
  #final_Omega_sparsity = prop_BSSC_nonlocal_Omega_1
  
  result = array(0, dim = c(p,p,D))
  
  for(i in 1:D){
    for(j in 1:p){
      for (k in 1:p){
        if (final_Omega_sparsity[,,i][j,k] == 0){
          result[,,i][j,k] = result[,,i][k,j] = 0
        } else {
          est_num = apply(Omega_iter_result[[i]]$Omega.mat[,,-(1:nburn)], c(1,2), sum, na.rm = TRUE)
          est_den = apply(Omega_iter_sparsity[[i]][,,(1:niter)], c(1,2), sum, na.rm = TRUE)
          result[,,i][j,k] = result[,,i][k,j] = est_num[j,k]/est_den[j,k]
        }
      }
      
    }
  }
  result = apply(result, c(1,2), mean, na.rm = TRUE)
  return(result)
}
