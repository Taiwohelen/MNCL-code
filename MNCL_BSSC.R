
#work directory and source file
setwd("")


#source("functions.R")
source("functions_github.R")

#libraries
library(truncnorm) #to generate truncnorn for refitted estimate
library(glasso)
library(MASS) # for mvrnorm function
library(matrixcalc) #to check positive definiteness
library(MLmetrics)
library(mltools)


####################################################
# Data generation
####################################################

p = 150
#p = 300
#p = 500
n = ceiling(p/2)
#n = ceiling(3*p/4)
mu = rep(0,p)


####################################################
# Generation of true precision matrix for X 
## The sparsity density is on the off diagonals(upper and lower combined)

Omega_start <- sparsity_generator(p, prop = 0.03, a = -0.6, b = -0.4, c = 0.4, d = 0.6)

is.positive.definite(Omega_start)
#Omega_start = Omega_start + 0.5*(diag(p))
is.positive.definite(Omega_start)

Omega0 = Omega_start

(length(which(Omega0!=0)) - p)/(p^2)   ##density of off diagonals

###################################################

# Generate D data sets of size n from MVN(0, inv(Omega0))
# X : n X p covariate matrix      ##this is the matrix of n p-variate observations

D <- 25    #no of dataset
X <- array(0, dim = c(n,p,D))
for(i in 1:D){
  X[,,i] <- mvrnorm(n, mu, Sigma = solve(Omega0))  ##works with the Sigma that is positive definite
}

####################################################

r = 1e-4    #hyperparameter for precision matrix
s = 1e-8     #hyperparameter for precision matrix

niter = 2000
nburn = 2000


############## MNCL Method

qn = 0.1
BSSC_nonlocal_Mode <- vector("list", D)
for(i in 1:D){
  BSSC_nonlocal_Mode[[i]] = BSSC_pMoM_Mode_max(X[,,i], qn, r, s, niter, nburn)
}

BSSC_nonlocal_Mode_time = rep(0,D)
for (i in 1:D){
  BSSC_nonlocal_Mode_time[[i]] <- mean(BSSC_nonlocal_Mode[[i]]$Time[-(1:nburn)]) #mean over niter for each dataset
}
BSSC_nonlocal_Mode_time_ave = mean(BSSC_nonlocal_Mode_time)

#########################

Omega0_classified = Omega0
Omega0_classified[which(Omega0_classified!=0)] = 1
Omega0_classified[which(Omega0_classified==0)] = 0

SP_BSSC_nonlocal_1 <- rep(0, D)
SE_BSSC_nonlocal_1 <- rep(0, D)
MCC_BSSC_nonlocal_1 <- rep(0, D)
BSSC_nonlocal_Omega_1 <- vector("list", D)
prop_BSSC_nonlocal_Omega_1 <- array(0, dim = c(p,p,D))
for (i in 1:D){
  
  #collecting the predicted precision matrix over the niter and no of replication
  BSSC_nonlocal_Omega_1[[i]] <- BSSC_nonlocal_Mode[[i]]$Omega.mat[,,-(1:nburn)]  ## direct Omega generated after burn-in for each iteration
  BSSC_nonlocal_Omega_1[[i]] <-  BSSC_nonlocal_Omega_1[[i]]!=0  #checks which entry is !=0 over all the iter
  prop_BSSC_nonlocal_Omega_1[,,i] <- apply(BSSC_nonlocal_Omega_1[[i]][,,(1:niter)], c(1, 2), mean, na.rm = TRUE)  #finds the mean over all the iter, returns one matrix
  
  #checking the sparsity pattern
  prop_BSSC_nonlocal_Omega_1[,,i][which(prop_BSSC_nonlocal_Omega_1[,,i] >= 0.5)] = 1
  prop_BSSC_nonlocal_Omega_1[,,i][which(prop_BSSC_nonlocal_Omega_1[,,i] < 0.5)] = 0

  # Accuracy metrics 
  SP_BSSC_nonlocal_1[i] = Specificity(Omega0_classified, prop_BSSC_nonlocal_Omega_1[,,i], positive = "1")
  SE_BSSC_nonlocal_1[i] = Sensitivity(Omega0_classified, prop_BSSC_nonlocal_Omega_1[,,i], positive = "1")
  MCC_BSSC_nonlocal_1[i] = ModelMetrics::mcc(Omega0_classified, prop_BSSC_nonlocal_Omega_1[,,i], cutoff = 0.5)
}

#mean value over the no of replications
SP_BSSC_nonlocal_1_ave = mean(SP_BSSC_nonlocal_1)
SE_BSSC_nonlocal_1_ave = mean(SE_BSSC_nonlocal_1)
MCC_BSSC_nonlocal_1_ave = mean(MCC_BSSC_nonlocal_1)

#estimate of the precision matrix over the number of replications
BSSC_pmom_estimate = BSSC_estimate(p, D, BSSC_nonlocal_Mode, BSSC_nonlocal_Omega_1, prop_BSSC_nonlocal_Omega_1)
BSSC_pmom_RE = frobenius.norm(BSSC_pmom_estimate - Omega0) / frobenius.norm(Omega0)


############ BSSC method
qn = 1/p
BSSC_norm <- vector("list", D)
for(i in 1:D){
  BSSC_norm[[i]] = BSSC_normal(X[,,i], qn, r, s, niter, nburn)
}


BSSC_norm_time <- rep(0, D)
for(i in 1:D){
  BSSC_norm_time[i] = mean(BSSC_norm[[i]]$Time[-(1:nburn)])   #mean over niter for each dataset
}
BSSC_normal_time_ave = mean(BSSC_norm_time)


#########


############# performance accuracy 

SP_BSSC_normal_1 <- rep(0, D)
SE_BSSC_normal_1 <- rep(0, D)
MCC_BSSC_normal_1 <- rep(0, D)
BSSC_normal_Omega_1 <- vector("list", D)
prop_BSSC_normal_Omega_1 <- array(0, dim = c(p,p,D))

for(i in 1:D){
  
  BSSC_normal_Omega_1[[i]] <- BSSC_norm[[i]]$Omega.mat[,,-(1:nburn)]
  BSSC_normal_Omega_1[[i]] <- BSSC_normal_Omega_1[[i]]!= 0
  prop_BSSC_normal_Omega_1[,,i] <- apply(BSSC_normal_Omega_1[[i]][,,(1:niter)], c(1,2), mean, na.rm = TRUE)
  
  prop_BSSC_normal_Omega_1[,,i][which(prop_BSSC_normal_Omega_1[,,i] >= 0.5)] = 1
  prop_BSSC_normal_Omega_1[,,i][which(prop_BSSC_normal_Omega_1[,,i] < 0.5)] = 0

  SP_BSSC_normal_1[i] = Specificity(Omega0_classified, prop_BSSC_normal_Omega_1[,,i], positive = "1")
  SE_BSSC_normal_1[i] = Sensitivity(Omega0_classified, prop_BSSC_normal_Omega_1[,,i], positive = "1")
  MCC_BSSC_normal_1[i] = ModelMetrics::mcc(Omega0_classified, prop_BSSC_normal_Omega_1[,,i], cutoff = 0.5)
  
}

SP_BSSC_normal_1_ave = mean(SP_BSSC_normal_1)
SE_BSSC_normal_1_ave = mean(SE_BSSC_normal_1)
MCC_BSSC_normal_1_ave = mean(MCC_BSSC_normal_1)

BSSC_normal_estimate = BSSC_estimate(p, D, BSSC_norm, BSSC_normal_Omega_1, prop_BSSC_normal_Omega_1)
BSSC_normal_RE = frobenius.norm(BSSC_normal_estimate - Omega0) / frobenius.norm(Omega0)

