# Generate four models

# Model1
Dynamic_Covariance2<- function(u,n){
  sigma_matrix<- matrix(0,n,n)
  u<-as.vector(as.matrix(u))
  for(i in 1:n){
    for(j in 1:n){
      sigma_matrix[i,j] <- exp(u[1])*(1/sqrt(2*pi)*exp(-(u[1])^2/2))^(abs(i-j))
    }
  }
  return(sigma_matrix)
}

# Model2
# Dynamic_Covariance2<- function(u,n){ 
#   sigma_matrix<- matrix(0,n,n)
#   u<-as.vector(as.matrix(u))
#   for(i in 1:n){
#     for(j in 1:n){
#       sigma_matrix[i,j] <- exp(u[1]+u[2])*(1/sqrt(2*pi)*exp(-((u[1]+u[2])/2)^2/2))^(abs(i-j))
#     }
#   }
#   return(sigma_matrix)
# }

# Model3
# Dynamic_Covariance2<- function(u,n){ 
#   sigma_matrix<- matrix(0,n,n)
#   u<-as.vector(as.matrix(u))
#   for(i in 1:n){
#     for(j in 1:n){
#       if(i==j){sigma_matrix[i,j]<-exp(2*u[1])}
#       else if((abs(i-j)==1) & (u[1]<=1 & u[1]>=-0.5)){sigma_matrix[i,j]<-exp(2*u[1])*0.5*exp(-(u[1]-0.25)^2/(0.75^2-(u[1]-0.25)^2))}
#       else if((abs(i-j)==2) & (u[1]<=1 & u[1]>= 0.3)){sigma_matrix[i,j]<-exp(2*u[1])*0.4*exp(-(u[1]-0.65)^2/(0.35^2-(u[1]-0.65)^2))}
#     }
#   }
#   return(sigma_matrix)
# }

# Model4
# Dynamic_Covariance<- function(u,n){ 
#   sigma_matrix<- matrix(0,n,n)
#   u<-as.vector(as.matrix(u))
#   for(i in 1:n){
#     for(j in 1:n){
#       if(i==j){sigma_matrix[i,j]<-exp(2*u[1])}
#       else if((abs(i-j)==1) & ((u[1]<=1 & u[1]>=-0.5)&(u[2]<=1 & u[2]>=-0.5))){sigma_matrix[i,j]<-exp(u[1]*2)*0.5*exp(-(u[1]-0.25)^2/(0.75^2-(u[1]-0.25)^2))}
#       else if((abs(i-j)==2) & ((u[1]<=1 & u[1]>= 0.3)&(u[2]<=1 & u[2]>= 0.3))){sigma_matrix[i,j]<-exp(u[1]*2)*0.4*exp(-(u[1]-0.65)^2/(0.35^2-(u[1]-0.65)^2))}
#     }
#   }
#   return(sigma_matrix)
# }
# Dynamic_Covariance2<- function(u,n){
#   u<-as.vector(as.matrix(u))
#   return(0.5*Dynamic_Covariance(c(u[1],u[2]),n)+0.5*Dynamic_Covariance(c(u[2],u[1]),n))
# }

# Mean-value function used for model1-4
Dynamic_mean<- function(u,n){
  mean_v<-rep(0,n)
  u<-as.vector(as.matrix(u))
  for(i in 1:n){
    mean_v[i]=0
  }
  return(mean_v)
}    

# SCAD threshold function
scad_penalty <- function(x, lambda, a) {
  abs_x <- abs(x)
  ifelse(abs_x <= 2*lambda, sign(x) *(abs_x-lambda)*((abs_x-lambda)>0), 
         ifelse(abs_x <= (a * lambda), ((a - 1) * x -sign(x)*a* lambda)/ (a-2), x))
}

# Kernel function
guassionkernel = function(x) exp(-0.5 * (x^2))/ sqrt(2 * pi)

# Main function
dcm = function(shushu, testpoint, X, Y, threshold.min = 0, threshold.max = 2) {
  # Sample size
  n  = length(X)
  # Dimension of Y
  pp = length(Y[1, ])
  error.NW = Y
  sigma.NW.jihe = list()
  
  ######################################
  #### Tuning the bandwidth parameter
  ######################################
  h1 = n^(-0.2) * median(abs(X - median(X)))/0.6745
  bandwidth = seq(h1, h1 + 0.9, by = 0.1) 
  
  NN = 10
  pppp = round(pp/12) 
  bandwidthA = matrix(0, nrow = NN, ncol = pppp)
  
  for (j in 1:NN) {
    a = sample(1:pp, pppp, replace = FALSE, prob = NULL)
    a = sort(a)
    bandwidthA[j, ] = a
  }
  
  error.band = numeric(length(bandwidth))
  
  for (band in 1:length(bandwidth)) {
    
    h = bandwidth[band]
    sum.NW = 0
    
    for (j in 1:NN) {
      
      a = bandwidthA[j, ]  
      
      YY = matrix(1, nrow = n, ncol = length(a))
      
      for (i in 1:length(a)) YY[, i]=Y[, a[i]]
      
      for (shu in 1:n) {
        x = X[shu]
        K = guassionkernel((X - x)/h)/h
        sigma.NW.jihe[[shu]] = matrix(0, ncol = length(a), nrow = length(a))
        for (i in 1:n) {
          if (i != shu) {
            sigma.NW.jihe[[shu]] = sigma.NW.jihe[[shu]] + K[i] * (t(t(YY[i, ])) %*% t(YY[i, ]))
          }
        }
        sigma.NW.jihe[[shu]] = sigma.NW.jihe[[shu]]/sum(K[-shu])
        sum.NW = sum.NW + t(YY[shu, ]) %*% solve(sigma.NW.jihe[[shu]]) %*% t(t(YY[shu, ])) + log(det(sigma.NW.jihe[[shu]]))
      }
      
    }
    error.band[band] = sum.NW/n/pp
  }
  
  h.wanzheng = bandwidth[which.min(error.band)]
  #############################################
  n.1 = round(0.7*n)
  n.2 = n - n.1
  
  CCC = h.wanzheng/n^(-0.2)
  h.train = CCC * n.1^(-0.2)
  h.test = CCC * n.2^(-0.2)
  
  # 30-fold for tuning dynamic thresholding parameter 
  N1 = 30
  TAa = matrix(0, nrow = N1, ncol = n.1)
  TAa.res = matrix(0, nrow = N1, ncol = n.2)
  
  for (j in 1:N1) {
    a = sample(1:n, n.1, replace = FALSE, prob = NULL)
    a = sort(a)
    res = (1:n)[is.na(pmatch(1:n, a))]
    TAa[j, ] = a
    TAa.res[j, ] = res
  }
  
  ##################################################################################################################
  ### Calculate the estimated covariance matrices for each test point.
  ##################################################################################################################
    
  x = testpoint[shushu]
  thresholding = seq(threshold.min, threshold.max, by = 0.1)
  
  K = guassionkernel((X - x)/h.wanzheng)/h.wanzheng
  sigma.NW.jihe = matrix(0, ncol = pp, nrow = pp)
  
  for (i in 1:n) {
    sigma.NW.jihe = sigma.NW.jihe + K[i] * (t(t(error.NW[i, ])) %*% t(error.NW[i, ]))
  }
  
  sigma.NW.jihe = sigma.NW.jihe/sum(K)
  sigma.z = numeric(pp)
  for (i in 1:n) {
    sigma.z = sigma.z + K[i] * error.NW[i, ]/sum(K)
  }
  
  sigma.NW.jihe = sigma.NW.jihe - t(t(sigma.z)) %*% t(sigma.z)
  diag_x=diag(sigma.NW.jihe)
  
  # The parameter used to do positive definiteness correction (Modified-Kernel)
  c0=0.001
  
  error.thresh.soft = numeric(length(thresholding))
  error.thresh.hard = numeric(length(thresholding))
  error.thresh.scad = numeric(length(thresholding))
  error.thresh.lasso = numeric(length(thresholding))
    
  for (thresh in 1:length(thresholding)) {
      
    s = thresholding[thresh]
    sum.NWs.soft = 0
    sum.NWs.hard = 0
    sum.NWs.scad = 0
    sum.NWs.lasso = 0
      
    for (j in 1:N1) {
        
      a = TAa[j, ]
      res = TAa.res[j, ]
      YYY = matrix(1, nrow = n.1, ncol = pp)
      YYY = Y[a, ]
      XXX = numeric(n.1)
      XXX = X[a]
      YYY.res = matrix(1, nrow = n.2, ncol = pp)
      YYY.res = Y[res, ]
      XXX.res = numeric(n.2)
      XXX.res = X[res]
        
      K = guassionkernel((XXX - x)/h.train)/h.train
        
      sigma.NW.jihe1 = matrix(0, ncol = pp, nrow = pp)
      sigma.NW.jihe.res1 = matrix(0, ncol = pp, nrow = pp)
        
      for (i in 1:n.1) {
        sigma.NW.jihe1 = sigma.NW.jihe1 + K[i] * (t(t(YYY[i, ])) %*% t(YYY[i, ]))
      }
      sigma.NW.jihe1 = sigma.NW.jihe1/sum(K)
      
      sigma.NW.jihe.soft1 = sign(sigma.NW.jihe1) *
        (abs(sigma.NW.jihe1) - s) *
        ((abs(sigma.NW.jihe1) - s) > 0)
      diag(sigma.NW.jihe.soft1)=diag_x
      
      sigma.NW.jihe.hard1 = (sigma.NW.jihe1) *((abs(sigma.NW.jihe1) - s) > 0)
      diag(sigma.NW.jihe.hard1)=diag_x
      
      sigma.NW.jihe.scad1 = scad_penalty(sigma.NW.jihe1,s,3.7)
      diag(sigma.NW.jihe.scad1)=diag_x
      
      sigma.NW.jihe.lasso1 = sign(sigma.NW.jihe1) *(abs(sigma.NW.jihe1) - s^2/abs(sigma.NW.jihe1))*((abs(sigma.NW.jihe1) - s^2/abs(sigma.NW.jihe1))  > 0)
      diag(sigma.NW.jihe.lasso1)=diag_x
      
      
      K.res = guassionkernel((XXX.res - x)/h.test)/h.test
        
      for (i in 1:n.2) {
        sigma.NW.jihe.res1 = sigma.NW.jihe.res1 + K.res[i] * (t(t(YYY.res[i,])) %*% t(YYY.res[i, ]))
      }
        
      sigma.NW.jihe.res1 = sigma.NW.jihe.res1/sum(K.res)
      
      HHH = sigma.NW.jihe.soft1 - sigma.NW.jihe.res1
      sum.NWs.soft = sum.NWs.soft + sum(diag(HHH %*% HHH))
      
      HHH = sigma.NW.jihe.hard1 - sigma.NW.jihe.res1
      sum.NWs.hard = sum.NWs.hard + sum(diag(HHH %*% HHH))
      
      HHH = sigma.NW.jihe.scad1 - sigma.NW.jihe.res1
      sum.NWs.scad = sum.NWs.scad + sum(diag(HHH %*% HHH))
      
      HHH = sigma.NW.jihe.lasso1 - sigma.NW.jihe.res1
      sum.NWs.lasso = sum.NWs.lasso + sum(diag(HHH %*% HHH))
        
    }
      
    error.thresh.soft[thresh] = sum.NWs.soft/n
    error.thresh.hard[thresh] = sum.NWs.hard/n
    error.thresh.scad[thresh] = sum.NWs.scad/n
    error.thresh.lasso[thresh] = sum.NWs.lasso/n
      
  }
    
  s.soft = thresholding[which.min(error.thresh.soft)]
  
  s.hard = thresholding[which.min(error.thresh.hard)]
  
  s.scad = thresholding[which.min(error.thresh.scad)]
  
  s.lasso = thresholding[which.min(error.thresh.lasso)]

  sigma.NW.jihe.soft = sign(sigma.NW.jihe) * (abs(sigma.NW.jihe) - s.soft) * ((abs(sigma.NW.jihe) - s.soft) > 0)
  diag(sigma.NW.jihe.soft)=diag_x
  a=min(eigen(sigma.NW.jihe.soft)$values)
  sigma.NW.jihe.soft_c=sigma.NW.jihe.soft+diag(rep(-a,pp)+c0)*(a<0)
  
  sigma.NW.jihe.hard = sigma.NW.jihe *((abs(sigma.NW.jihe) - s.hard) > 0)
  diag(sigma.NW.jihe.hard)=diag_x
  a=min(eigen(sigma.NW.jihe.hard)$values)
  sigma.NW.jihe.hard_c=sigma.NW.jihe.hard+diag(rep(-a,pp)+c0)*(a<0)
  
  sigma.NW.jihe.scad = scad_penalty(sigma.NW.jihe,s.scad,3.7)
  diag(sigma.NW.jihe.scad)=diag_x
  a=min(eigen(sigma.NW.jihe.scad)$values)
  sigma.NW.jihe.scad_c=sigma.NW.jihe.scad+diag(rep(-a,pp)+c0)*(a<0)
  
  sigma.NW.jihe.lasso = sign(sigma.NW.jihe) *(abs(sigma.NW.jihe) - s.lasso^2/abs(sigma.NW.jihe))*((abs(sigma.NW.jihe) - s.lasso^2/abs(sigma.NW.jihe)) > 0)
  diag(sigma.NW.jihe.lasso)=diag_x
  a=min(eigen(sigma.NW.jihe.lasso)$values)
  sigma.NW.jihe.lasso_c=sigma.NW.jihe.lasso+diag(rep(-a,pp)+c0)*(a<0)
    
  
  return(list(sigma.NW.jihe.soft,sigma.NW.jihe.hard,sigma.NW.jihe.scad,sigma.NW.jihe.lasso,sigma.NW.jihe.soft_c,sigma.NW.jihe.hard_c,sigma.NW.jihe.scad_c,sigma.NW.jihe.lasso_c))
  
}



library(MASS)
library(lava)
library(BiocParallel)
library(parallel)

# Number of cores
ncores = 40 
mcparam = SnowParam(workers = ncores)
# Number of repetitions
NN=50 

# MFL,MSL for Kernel
MFL_soft=matrix(0, ncol = NN, nrow = 1)
MSL_soft=matrix(0, ncol = NN, nrow = 1)

MFL_hard=matrix(0, ncol = NN, nrow = 1)
MSL_hard=matrix(0, ncol = NN, nrow = 1)

MFL_lasso=matrix(0, ncol = NN, nrow = 1)
MSL_lasso=matrix(0, ncol = NN, nrow = 1)

MFL_scad=matrix(0, ncol = NN, nrow = 1)
MSL_scad=matrix(0, ncol = NN, nrow = 1)


# MFL,MSL for Modified-Kernel 
MFL_soft_c=matrix(0, ncol = NN, nrow = 1)
MSL_soft_c=matrix(0, ncol = NN, nrow = 1)

MFL_hard_c=matrix(0, ncol = NN, nrow = 1)
MSL_hard_c=matrix(0, ncol = NN, nrow = 1)

MFL_lasso_c=matrix(0, ncol = NN, nrow = 1)
MSL_lasso_c=matrix(0, ncol = NN, nrow = 1)

MFL_scad_c=matrix(0, ncol = NN, nrow = 1)
MSL_scad_c=matrix(0, ncol = NN, nrow = 1)


# TPR,FPR for Kernel (only used for model 3 and 4)
TPR_soft=matrix(0, ncol = NN, nrow = 1)
FPR_soft=matrix(0, ncol = NN, nrow = 1)

TPR_hard=matrix(0, ncol = NN, nrow = 1)
FPR_hard=matrix(0, ncol = NN, nrow = 1)

TPR_lasso=matrix(0, ncol = NN, nrow = 1)
FPR_lasso=matrix(0, ncol = NN, nrow = 1)

TPR_scad=matrix(0, ncol = NN, nrow = 1)
FPR_scad=matrix(0, ncol = NN, nrow = 1)



for (ii in 1:NN){
  # Sample size
  n <- 100 
  # Dimension of U
  p1 <- 10 #20
  # Dimension of Y
  p2 <- 100 #150,200
  
  # Generate data
  set.seed(ii)
  X <- matrix(runif(n * p1)*2-1, n, p1)
  Y <- matrix(0,n,p2)
  for (i in 1:n)
  {
    set.seed(i)
    Y[i,]<-mvrnorm(1, Dynamic_mean(X[i,],p2), Dynamic_Covariance2(X[i,],p2))
  }
  
  # Generate 30 test points
  ns=30
  set.seed(ii+1) # used for model 1,4; set.seed(ii+6) for model 2,3
  testpoint<- matrix(runif(ns * p1)*2-1, ns,p1)
  ##########################################################################################################
  ##########################################################################################################
  
  thresh1 = numeric(ns)
  
  # Result for covariate u_1 (used for model 1-4)
  results = bplapply(seq(1:ns), dcm, testpoint[,1], X[,1], Y ,BPPARAM = mcparam)
  # Result for covariate u_2 (used for model 1,3)
  # results = bplapply(seq(1:ns), dcm, testpoint[,2], X[,2], Y ,BPPARAM = mcparam)
  # Result for covariate u_3 (used for model 2,4)
  # results = bplapply(seq(1:ns), dcm, testpoint[,3], X[,3], Y ,BPPARAM = mcparam)
  
  F_norm_soft<-matrix(0, ncol =ns, nrow = 1)
  F_norm_hard<-matrix(0, ncol =ns, nrow = 1)
  F_norm_scad<-matrix(0, ncol =ns, nrow = 1)
  F_norm_lasso<-matrix(0, ncol =ns, nrow = 1)
  F_norm_soft_1<-matrix(0, ncol =ns, nrow = 1)
  F_norm_hard_1<-matrix(0, ncol =ns, nrow = 1)
  F_norm_scad_1<-matrix(0, ncol =ns, nrow = 1)
  F_norm_lasso_1<-matrix(0, ncol =ns, nrow = 1)
  
  spectral_radius_soft<-matrix(0, ncol =ns, nrow = 1)
  spectral_radius_hard<-matrix(0, ncol =ns, nrow = 1)
  spectral_radius_scad<-matrix(0, ncol =ns, nrow = 1)
  spectral_radius_lasso<-matrix(0, ncol =ns, nrow = 1)
  spectral_radius_soft_1<-matrix(0, ncol =ns, nrow = 1)
  spectral_radius_hard_1<-matrix(0, ncol =ns, nrow = 1)
  spectral_radius_scad_1<-matrix(0, ncol =ns, nrow = 1)
  spectral_radius_lasso_1<-matrix(0, ncol =ns, nrow = 1)
  
  # Only used for model 3 and 4
  FPR_loss_soft<-matrix(0, ncol =ns, nrow = 1)
  FPR_loss_hard<-matrix(0, ncol =ns, nrow = 1)
  FPR_loss_scad<-matrix(0, ncol =ns, nrow = 1)
  FPR_loss_lasso<-matrix(0, ncol =ns, nrow = 1)
  TPR_loss_soft<-matrix(0, ncol =ns, nrow = 1)
  TPR_loss_hard<-matrix(0, ncol =ns, nrow = 1)
  TPR_loss_scad<-matrix(0, ncol =ns, nrow = 1)
  TPR_loss_lasso<-matrix(0, ncol =ns, nrow = 1)
  
  for(kk in 1:ns)
  {
    eig <- eigen(results[[kk]][[1]]-Dynamic_Covariance2(testpoint[kk,],p2))
    eig_vals <- eig$values
    spectral_radius_soft[1,kk] <- max(abs(eig_vals))
    F_norm_soft[1,kk]<-sqrt(tr((results[[kk]][[1]]-Dynamic_Covariance2(testpoint[kk,],p2))
                               %*%(results[[kk]][[1]]-Dynamic_Covariance2(testpoint[kk,],p2))))
    
    
    eig <- eigen(results[[kk]][[2]]-Dynamic_Covariance2(testpoint[kk,],p2))
    eig_vals <- eig$values
    spectral_radius_hard[1,kk] <- max(abs(eig_vals))
    F_norm_hard[1,kk]<-sqrt(tr((results[[kk]][[2]]-Dynamic_Covariance2(testpoint[kk,],p2))
                               %*%(results[[kk]][[2]]-Dynamic_Covariance2(testpoint[kk,],p2))))
    
    eig <- eigen(results[[kk]][[3]]-Dynamic_Covariance2(testpoint[kk,],p2))
    eig_vals <- eig$values
    spectral_radius_scad[1,kk] <- max(abs(eig_vals))
    F_norm_scad[1,kk]<-sqrt(tr((results[[kk]][[3]]-Dynamic_Covariance2(testpoint[kk,],p2))
                               %*%(results[[kk]][[3]]-Dynamic_Covariance2(testpoint[kk,],p2))))
    
    eig <- eigen(results[[kk]][[4]]-Dynamic_Covariance2(testpoint[kk,],p2))
    eig_vals <- eig$values
    spectral_radius_lasso[1,kk] <- max(abs(eig_vals))
    F_norm_lasso[1,kk]<-sqrt(tr((results[[kk]][[4]]-Dynamic_Covariance2(testpoint[kk,],p2))
                                %*%(results[[kk]][[4]]-Dynamic_Covariance2(testpoint[kk,],p2))))
    
    
    
    eig <- eigen(results[[kk]][[5]]-Dynamic_Covariance2(testpoint[kk,],p2))
    eig_vals <- eig$values
    spectral_radius_soft_1[1,kk] <- max(abs(eig_vals))
    F_norm_soft_1[1,kk]<-sqrt(tr((results[[kk]][[5]]-Dynamic_Covariance2(testpoint[kk,],p2))
                                 %*%(results[[kk]][[5]]-Dynamic_Covariance2(testpoint[kk,],p2))))
    
    
    eig <- eigen(results[[kk]][[6]]-Dynamic_Covariance2(testpoint[kk,],p2))
    eig_vals <- eig$values
    spectral_radius_hard_1[1,kk] <- max(abs(eig_vals))
    F_norm_hard_1[1,kk]<-sqrt(tr((results[[kk]][[6]]-Dynamic_Covariance2(testpoint[kk,],p2))
                                 %*%(results[[kk]][[6]]-Dynamic_Covariance2(testpoint[kk,],p2))))
    
    eig <- eigen(results[[kk]][[7]]-Dynamic_Covariance2(testpoint[kk,],p2))
    eig_vals <- eig$values
    spectral_radius_scad_1[1,kk] <- max(abs(eig_vals))
    F_norm_scad_1[1,kk]<-sqrt(tr((results[[kk]][[7]]-Dynamic_Covariance2(testpoint[kk,],p2))
                                 %*%(results[[kk]][[7]]-Dynamic_Covariance2(testpoint[kk,],p2))))
    
    eig <- eigen(results[[kk]][[8]]-Dynamic_Covariance2(testpoint[kk,],p2))
    eig_vals <- eig$values
    spectral_radius_lasso_1[1,kk] <- max(abs(eig_vals))
    F_norm_lasso_1[1,kk]<-sqrt(tr((results[[kk]][[8]]-Dynamic_Covariance2(testpoint[kk,],p2))
                                  %*%(results[[kk]][[8]]-Dynamic_Covariance2(testpoint[kk,],p2))))
    
    FPR_loss_soft[1,kk]<-sum(results[[kk]][[1]]!=0 & Dynamic_Covariance2(testpoint[kk,],p2) ==0)/sum(Dynamic_Covariance2(testpoint[kk,],p2) ==0)
    FPR_loss_hard[1,kk]<-sum(results[[kk]][[2]]!=0 & Dynamic_Covariance2(testpoint[kk,],p2) ==0)/sum(Dynamic_Covariance2(testpoint[kk,],p2) ==0)
    FPR_loss_scad[1,kk]<-sum(results[[kk]][[3]]!=0 & Dynamic_Covariance2(testpoint[kk,],p2) ==0)/sum(Dynamic_Covariance2(testpoint[kk,],p2) ==0)
    FPR_loss_lasso[1,kk]<-sum(results[[kk]][[4]]!=0 & Dynamic_Covariance2(testpoint[kk,],p2) ==0)/sum(Dynamic_Covariance2(testpoint[kk,],p2) ==0)
    
    TPR_loss_soft[1,kk]<-sum(results[[kk]][[1]]!=0 & Dynamic_Covariance2(testpoint[kk,],p2) !=0)/sum(Dynamic_Covariance2(testpoint[kk,],p2) !=0)
    TPR_loss_hard[1,kk]<-sum(results[[kk]][[2]]!=0 & Dynamic_Covariance2(testpoint[kk,],p2) !=0)/sum(Dynamic_Covariance2(testpoint[kk,],p2) !=0)
    TPR_loss_scad[1,kk]<-sum(results[[kk]][[3]]!=0 & Dynamic_Covariance2(testpoint[kk,],p2) !=0)/sum(Dynamic_Covariance2(testpoint[kk,],p2) !=0)
    TPR_loss_lasso[1,kk]<-sum(results[[kk]][[4]]!=0 & Dynamic_Covariance2(testpoint[kk,],p2) !=0)/sum(Dynamic_Covariance2(testpoint[kk,],p2) !=0)
  }
  MFL_soft[ii]=median(F_norm_soft)
  MSL_soft[ii]=median(spectral_radius_soft)
  
  MFL_hard[ii]=median(F_norm_hard)
  MSL_hard[ii]=median(spectral_radius_hard)
  
  MFL_scad[ii]=median(F_norm_scad)
  MSL_scad[ii]=median(spectral_radius_scad)
  
  MFL_lasso[ii]=median(F_norm_lasso)
  MSL_lasso[ii]=median(spectral_radius_lasso)
  
  
  MFL_soft_c[ii]=median(F_norm_soft_1)
  MSL_soft_c[ii]=median(spectral_radius_soft_1)
  
  MFL_hard_c[ii]=median(F_norm_hard_1)
  MSL_hard_c[ii]=median(spectral_radius_hard_1)
  
  MFL_scad_c[ii]=median(F_norm_scad_1)
  MSL_scad_c[ii]=median(spectral_radius_scad_1)
  
  MFL_lasso_c[ii]=median(F_norm_lasso_1)
  MSL_lasso_c[ii]=median(spectral_radius_lasso_1)
  
  TPR_soft[ii]=median(TPR_loss_soft)
  FPR_soft[ii]=median(FPR_loss_soft)
  
  TPR_hard[ii]=median(TPR_loss_hard)
  FPR_hard[ii]=median(FPR_loss_hard)
  
  TPR_scad[ii]=median(TPR_loss_scad)
  FPR_scad[ii]=median(FPR_loss_scad)
  
  TPR_lasso[ii]=median(TPR_loss_lasso)
  FPR_lasso[ii]=median(FPR_loss_lasso)
}

# Compute mean of MFL,MSL,TPR and FPR under 50 repetitions
MFL_soft_mean <- mean(MFL_soft)
MSL_soft_mean <- mean(MSL_soft)
MFL_hard_mean <- mean(MFL_hard)
MSL_hard_mean <- mean(MSL_hard)
MFL_lasso_mean <- mean(MFL_lasso)
MSL_lasso_mean <- mean(MSL_lasso)
MFL_scad_mean <- mean(MFL_scad)
MSL_scad_mean <- mean(MSL_scad)

MFL_soft_c_mean <- mean(MFL_soft_c)
MSL_soft_c_mean <- mean(MSL_soft_c)
MFL_hard_c_mean <- mean(MFL_hard_c)
MSL_hard_c_mean <- mean(MSL_hard_c)
MFL_lasso_c_mean <- mean(MFL_lasso_c)
MSL_lasso_c_mean <- mean(MSL_lasso_c)
MFL_scad_c_mean <- mean(MFL_scad_c)
MSL_scad_c_mean <- mean(MSL_scad_c)

TPR_soft_mean <- mean(TPR_soft)
FPR_soft_mean <- mean(FPR_soft)
TPR_hard_mean <- mean(TPR_hard)
FPR_hard_mean <- mean(FPR_hard)
TPR_lasso_mean <- mean(TPR_lasso)
FPR_lasso_mean <- mean(FPR_lasso)
TPR_scad_mean <- mean(TPR_scad)
FPR_scad_mean <- mean(FPR_scad)

# Compute sd of MFL,MSL,TPR and FPR under 50 repetitions
MFL_soft_sd <- sd(MFL_soft)
MSL_soft_sd <- sd(MSL_soft)
MFL_hard_sd <- sd(MFL_hard)
MSL_hard_sd <- sd(MSL_hard)
MFL_lasso_sd <- sd(MFL_lasso)
MSL_lasso_sd <- sd(MSL_lasso)
MFL_scad_sd <- sd(MFL_scad)
MSL_scad_sd <- sd(MSL_scad)

MFL_soft_c_sd <- sd(MFL_soft_c)
MSL_soft_c_sd <- sd(MSL_soft_c)
MFL_hard_c_sd <- sd(MFL_hard_c)
MSL_hard_c_sd <- sd(MSL_hard_c)
MFL_lasso_c_sd <- sd(MFL_lasso_c)
MSL_lasso_c_sd <- sd(MSL_lasso_c)
MFL_scad_c_sd <- sd(MFL_scad_c)
MSL_scad_c_sd <- sd(MSL_scad_c)

TPR_soft_sd <- sd(TPR_soft)
FPR_soft_sd <- sd(FPR_soft)
TPR_hard_sd <- sd(TPR_hard)
FPR_hard_sd <- sd(FPR_hard)
TPR_lasso_sd <- sd(TPR_lasso)
FPR_lasso_sd <- sd(FPR_lasso)
TPR_scad_sd <- sd(TPR_scad)
FPR_scad_sd <- sd(FPR_scad)


x1 <- paste(sprintf("%0.2f",MFL_hard_mean),"(",sprintf("%0.2f",MFL_hard_sd),")","&",
            sprintf("%0.2f",MFL_lasso_mean),"(",sprintf("%0.2f",MFL_lasso_sd),")","&",
            sprintf("%0.2f",MFL_scad_mean),"(",sprintf("%0.2f",MFL_scad_sd),")","&&",
            sprintf("%0.2f",MFL_soft_mean),"(",sprintf("%0.2f",MFL_soft_sd),")",sep="")

x2 <- paste(sprintf("%0.2f",MSL_hard_mean),"(",sprintf("%0.2f",MSL_hard_sd),")","&",
            sprintf("%0.2f",MSL_lasso_mean),"(",sprintf("%0.2f",MSL_lasso_sd),")","&",
            sprintf("%0.2f",MSL_scad_mean),"(",sprintf("%0.2f",MSL_scad_sd),")","&&",
            sprintf("%0.2f",MSL_soft_mean),"(",sprintf("%0.2f",MSL_soft_sd),")",sep="")

x3 <- paste(sprintf("%0.2f",MFL_hard_c_mean),"(",sprintf("%0.2f",MFL_hard_c_sd),")","&",
            sprintf("%0.2f",MFL_lasso_c_mean),"(",sprintf("%0.2f",MFL_lasso_c_sd),")","&",
            sprintf("%0.2f",MFL_scad_c_mean),"(",sprintf("%0.2f",MFL_scad_c_sd),")","&&",
            sprintf("%0.2f",MFL_soft_c_mean),"(",sprintf("%0.2f",MFL_soft_c_sd),")",sep="")

x4 <- paste(sprintf("%0.2f",MSL_hard_c_mean),"(",sprintf("%0.2f",MSL_hard_c_sd),")","&",
            sprintf("%0.2f",MSL_lasso_c_mean),"(",sprintf("%0.2f",MSL_lasso_c_sd),")","&",
            sprintf("%0.2f",MSL_scad_c_mean),"(",sprintf("%0.2f",MSL_scad_c_sd),")","&&",
            sprintf("%0.2f",MSL_soft_c_mean),"(",sprintf("%0.2f",MSL_soft_c_sd),")",sep="")

x5 <- paste(sprintf("%0.2f",TPR_hard_mean),"(",sprintf("%0.2f",TPR_hard_sd),")","&",
            sprintf("%0.2f",TPR_lasso_mean),"(",sprintf("%0.2f",TPR_lasso_sd),")","&",
            sprintf("%0.2f",TPR_scad_mean),"(",sprintf("%0.2f",TPR_scad_sd),")","&&",
            sprintf("%0.2f",TPR_soft_mean),"(",sprintf("%0.2f",TPR_soft_sd),")",sep="")

x6 <- paste(sprintf("%0.2f",FPR_hard_mean),"(",sprintf("%0.2f",FPR_hard_sd),")","&",
            sprintf("%0.2f",FPR_lasso_mean),"(",sprintf("%0.2f",FPR_lasso_sd),")","&",
            sprintf("%0.2f",FPR_scad_mean),"(",sprintf("%0.2f",FPR_scad_sd),")","&&",
            sprintf("%0.2f",FPR_soft_mean),"(",sprintf("%0.2f",FPR_soft_sd),")",sep="")

# Print the result
x1
x2
x3
x4
x5
x6

