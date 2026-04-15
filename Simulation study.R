rm(list=ls(all=TRUE))
start_time <- Sys.time()

# Set seed for reproducibility
set.seed(123)
 
 
# Load required libraries

library(doParallel)
library(foreach)
 
# Initialize variables
p <-500
n1 <-50
n2 <-50

alpha <- 0.05

########################################## Covariance Structure ############################

#### AR(1), Non-sparse Dependence Structure , and LRD Dependence Structure

################## Start: AR (1)   ##################

#Sigma<- toeplitz(0.6^(0:(p-1)))

################## End: AR (1)   ##################



############### Start:Non-sparse Dependence Structure ############################

#Sigma<-matrix(NA,p,p)
#for (i in 1:p){
   # for (j in 1:p){
      #Sigma[i,j] = abs(i-j)^{-5}/2
    #}
  #}
  #for (i in 1:p){
    #Sigma[i,i] = 1
  #}

############### End :Non-sparse Dependence Structure ############################

############ Start: LRD structure  #############################

H <- 0.60  # Hurst parameter
Sigma <- matrix(0, p, p)
for (i in 1:p) {
  for (j in 1:p) {
    k1 <- abs(i - j)
    Sigma[i, j] <- 0.5 * ((k1 + 1)^(2*H) + (k1 - 1)^(2*H) - 2*k1^(2*H))
  }
}

diag(Sigma)<-1
############ Start: LRD structure  #############################


SS<-chol(Sigma)



############ parameter values for power 

a<-0.3  # signals sparsity parameter. It considered two values, 0.3 and 0.45 in the simulation

l<-1/2  # determines the strength of signals

r2<-ceiling(p^{1-a})  ####### no. of positive signals

delta1<-(log(p))^{1/2}/(n1+n2)^{l} # values of delta under H1
delta0<-0  ###### under null hypothesis

diff <- c()
for (i in 1:r2) {
  diff[i] <- i
}


#mu1 <-matrix(rep(c(rep(delta0,r2), rep(0,(p-r2)))),nrow=n1, ncol=p, byrow=T)
mu1 <-matrix(rep(c(rep(delta1,r2), rep(0,(p-r2)))),nrow=n1, ncol=p, byrow=T)  

#L<- ceiling(p^{1/2})  ### For Week Dependence
L <- ceiling(2*p^{1/2})  ### For Long Range Dependence

# Prepare for parallel processing
cl <- makeCluster(detectCores() - 5)  # Reserve one core for system stability
registerDoParallel(cl)
start <- Sys.time()
# Parallel execution using 'foreach'
results <- foreach(j = 1:1000, .combine = rbind) %dopar% {
  #library(expm)  # For matrix operations
  library(Matrix)
  #library(evd)
  library(goftest)
  #library(MASS)
  library(Matrix)
  library(mvtnorm)
  library(highmean)
  library(Rfast)
 
  library(Matrix)
  library(stats)
  library(matrixcalc)  
  # Here you will insert the computations that need to be parallelized.
  # This often includes data generation, application of statistical tests,
  # and collection of results. For example:
 
f<- function(u){
  return(u*(1/4 + 1/(2*pi)*asin(u)) + (sqrt(1 - u^2) - 1)/(2*pi))
}




 
  # Simulate data
   
 #X <- mu1 + t(SS %*% matrix(rgamma(n1*p,2,1)-2, nrow=p, ncol =n1))
 #Y <- t(SS %*% matrix(rgamma(n2*p,2,1)-2, nrow=p, ncol =n2))

 X <- mu1 + t(SS %*% matrix(rnorm(n1*p,0,1), nrow=p, ncol =n1))
  Y <- t(SS %*% matrix(rnorm(n2*p,0,1), nrow=p, ncol =n2))

   #X <- mu1 + t(SS %*% matrix(rt(n1 * p,5), ncol = n1))
   #Y <- t(SS %*% matrix(rt(n2 * p,5), ncol = n2))

 
  W <- mean(apply(data.frame(cbind((colMeans(X) - colMeans(Y))/sqrt(apply(X, 2, var)/n1+ apply(Y, 2, var)/n2), rep(0,p))), 1, max))
 
  ### mu Hat of W under the null hypothesis
  WMu<- 1/2 * sqrt(2/pi)
 
  ########################## vriance estimate of the W ##############################################
 
  eta.matrix<-diag(sqrt(1/diag(cov(X)/n1+cov(Y)/n2)),nrow=p)%*%(cov(X)/n1+cov(Y)/n2)%*%diag(sqrt(1/diag(cov(X)/n1+cov(Y)/n2)),nrow=p) ###### matrix with components rho^{*}_{ij}
  
  diag(eta.matrix)<-0
  
  gamma.hat<-apply(eta.matrix,2,f) ####### computes gamma.hat matrix with diagonal elemnst 0
  
  diag(gamma.hat)<-1/2*(1-1/pi)  ####### computes gamma.hat matrix 
 
 band_mask <- abs(row(gamma.hat) - col(gamma.hat)) <= L
 WVar <- sum(gamma.hat[band_mask]) / p
  
 #################### SMC Test #############################
 
 
Tobs<- sqrt(p)*(W - WMu) / sqrt(WVar)
 pvals <- 1 - pnorm(Tobs, 0, 1)
if (pvals <alpha) count1<-1 else count1<-0
 
 
################ Chongcharoen and others ##################################


 make_block_sizes <- function(p, v, q = NULL) {
  # Recommended in the paper: q = v - 6
  if (is.null(q)) q <- v - 6

  if (!is.numeric(q) || length(q) != 1 || q <= 0) {
    stop("q must be a positive integer.")
  }
  q <- as.integer(q)

  if (q >= v - 1) {
    stop("Need q < v - 1 so each block covariance is invertible.")
  }
  if (q >= v - 3) {
    warning("For the variance formula, it is safer to use q < v - 3.")
  }

  n_full <- p %/% q
  rem <- p %% q

  if (rem == 0) {
    qs <- rep(q, n_full)
  } else {
    qs <- c(rep(q, n_full), rem)
  }

  if (any(qs >= v - 1)) {
    stop("At least one block size is too large: need every q_j < v - 1.")
  }
  qs
}

# ---------- helper: pooled covariance ----------
pooled_cov <- function(X, Y) {
  n1 <- nrow(X)
  n2 <- nrow(Y)
  v  <- n1 + n2 - 2
  q <- v - 6

  S1 <- stats::cov(X)   # denominator n1 - 1
  S2 <- stats::cov(Y)   # denominator n2 - 1

  Sp <- ((n1 - 1) * S1 + (n2 - 1) * S2) / v
  list(Sp = Sp, S1 = S1, S2 = S2, v = v)
}

# ---------- main test function ----------

     v  <- n1 + n2 - 2
    q <- v - 6

  if (p <= 1) stop("Need p >= 2.")
  if (n1 < 2 || n2 < 2) stop("Need n1, n2 >= 2.")

  qs <- make_block_sizes(p = p, v = v, q = q)

  covs <- pooled_cov(X, Y)
  Sp   <- covs$Sp

  xbar1 <- colMeans(X)
  xbar2 <- colMeans(Y)
  dbar  <- xbar1 - xbar2

  # Build block-diagonal inverse via blockwise inversion
  Tn <- 0
  start <- 1

  for (qj in qs) {
    idx <- start:(start + qj - 1)
    Spj <- Sp[idx, idx, drop = FALSE]
    dbj <- dbar[idx]

    # Use solve; qr.solve could be used if numerical issues arise
    Tn <- Tn + as.numeric(t(dbj) %*% solve(Spj, dbj))

    start <- start + qj
  }

  Tn <- (n1 * n2 / (n1 + n2)) * Tn

  # Mean and variance correction from the paper
  mean_corr <- sum(v * qs / (v - qs - 1))
  var_corr  <- sum(
    2 * v^2 * (v - 1) * qs /
      ((v - qs + 1) * (v - qs - 1)^2 * (v - qs - 3))
  )

  Tq <- (Tn - mean_corr) / sqrt(var_corr)

  # One-sided Follmann-style sign condition
  sign_stat <- sum(dbar)

  crit <- stats::qnorm(1 - 2 * alpha)
  reject <- (Tq >= crit) && (sign_stat > 0)
    if((Tq >= crit) & (sign_stat >0)) count2<-1 else count2<-0
 



###########################################################
  SD<-apval_Sri2008(X, Y)
  CQ<-apval_Chen2010(X,Y)
  CXL<-apval_Cai2014(X,Y)
 
 
 
  if((SD$pval< 2*alpha) & (sum(colMeans(X) - colMeans(Y))>0)) count3<-1 else count3<-0
  if((CQ$pval< 2*alpha)& (sum(colMeans(X) - colMeans(Y))>0)) count4<-1 else count4<-0
  if((CXL$pval< 2*alpha) & (sum(colMeans(X) - colMeans(Y))>0)) count5<-1 else count5<-0
  return(c(count1, count2, count3, count4, count5))
  print(j)
}

# Stop the parallel cluster
stopCluster(cl)

# Calculate overall results
mean_results <- colMeans(results[,1:5])
print(mean_results)


print(paste("Empirical Type I Error rate of SMC:",mean_results[[1]]))
print(paste("Empirical Type I Error rate of SC:",mean_results[[2]]))
print(paste("Empirical Type I Error rate of SD+:",mean_results[[3]]))
print(paste("Empirical Type I Error rate of CQ+:",mean_results[[4]]))
print(paste("Empirical Type I Error rate of Cai+:",mean_results[[5]]))

##################
print(paste("Empirical power of SMC:",mean_results[[1]]))
print(paste("Empirical power of SC:",mean_results[[2]]))
print(paste("Empirical power of SD+:",mean_results[[3]]))
print(paste("Empirical power of CQ+:",mean_results[[4]]))
print(paste("Empirical power of Cai+:",mean_results[[5]]))

##### 

# Calculate elapsed time
end_time <- Sys.time()
print(end_time - start_time)






