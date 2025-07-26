# Generate four models

# Mode1 d=5
Dynamic_Covariance2<- function(u,n){
  sigma_matrix<- matrix(0,n,n)
  u<-as.vector(as.matrix(u))
  for(i in 1:n){
    for(j in 1:n){
      sigma_matrix[i,j] <- exp(sum(c(3,rep(0.1,3),0)*u))*(1/sqrt(2*pi)*exp(-(sum(c(30/33,rep(1/33,3),0)*u))^2/2))^(abs(i-j))
    }
  }
  return(sigma_matrix)
}

# Mode2 d=10
# Dynamic_Covariance2<- function(u,n){
#   sigma_matrix<- matrix(0,n,n)
#   u<-as.vector(as.matrix(u))
#   for(i in 1:n){
#     for(j in 1:n){
#       sigma_matrix[i,j] <- exp(sum(c(3,rep(0.1,6),rep(0,3))*u))*(1/sqrt(2*pi)*exp(-(sum(c(30/36,rep(1/36,6),rep(0,3))*u))^2/2))^(abs(i-j))
#     }
#   }
#   return(sigma_matrix)
# }

# Model3 d=5
# Dynamic_Covariance2<- function(u,n){
#   sigma_matrix<- matrix(0,n,n)
#   u<-as.vector(as.matrix(u))
#   for(i in 1:n){
#     for(j in 1:n){
#       if(i==j){sigma_matrix[i,j]<-exp(sum(c(3,rep(0.1,3),0)*u))}
#       else if((abs(i-j)==1) & ((u[1]<=0.9 & u[1]>=-0.5)&(u[2]<=0.9 & u[2]>=-0.5))){
#         sigma_matrix[i,j]<-exp(sum(c(3,rep(0.1,3),0)*u))*0.5*exp(-(u[1]-0.25)^2/(0.75^2-(u[1]-0.25)^2))}
#       else if((abs(i-j)==2) & ((u[1]<=0.9 & u[1]>= 0.3)&(u[2]<=0.9 & u[2]>= 0.3))){
#         sigma_matrix[i,j]<-exp(sum(c(3,rep(0.1,3),0)*u))*0.4*exp(-(u[1]-0.65)^2/(0.35^2-(u[1]-0.65)^2))}
#     }
#   }
#   return(sigma_matrix)
# }

# Model4 d=10
# Dynamic_Covariance2<- function(u,n){
#   sigma_matrix<- matrix(0,n,n)
#   u<-as.vector(as.matrix(u))
#   for(i in 1:n){
#     for(j in 1:n){
#       if(i==j){sigma_matrix[i,j]<-exp(sum(c(3,rep(0.1,6),rep(0,3))*u))}
#       else if((abs(i-j)==1) & ((u[1]<=0.9 & u[1]>=-0.5)&(u[2]<=0.9 & u[2]>=-0.5))){
#         sigma_matrix[i,j]<-exp(sum(c(3,rep(0.1,6),rep(0,3))*u))*0.5*exp(-(u[1]-0.25)^2/(0.75^2-(u[1]-0.25)^2))}
#       else if((abs(i-j)==2) & ((u[1]<=0.9 & u[1]>= 0.3)&(u[2]<=0.9 & u[2]>= 0.3))){
#         sigma_matrix[i,j]<-exp(sum(c(3,rep(0.1,6),rep(0,3))*u))*0.4*exp(-(u[1]-0.65)^2/(0.35^2-(u[1]-0.65)^2))}
#     }
#   }
#   return(sigma_matrix)
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

# Multivariate Kernel function
guassionkernel = function(x) {
  kernelmatrix <- exp(-0.5 * (x^2))/ sqrt(2 * pi)
  result <- apply(kernelmatrix,1,prod)
  return(result)
}

# Main function
dcm = function(shushu, testpoint, X, Y, threshold.min = 0, threshold.max = 2) {
  # Sample size
  n  = length(X[,1])
  # Dimension of X
  px = length(X[1,])
  # Dimension of Y
  pp = length(Y[1, ])
  error.NW = Y
  sigma.NW.jihe = list()
  
  ######################################
  #### Tunning the bandwidth parameter
  ######################################
  h.wanzheng = c(0.27,1.16,1.16,1.16,1.16)^(5/(4+px)) #d=5
  # h.wanzheng = c(0.27,1.16,1.16,rep(1.16,7))^(5/(4+px)) #d=10
  #############################################
  n.1 = round(0.7*n)
  n.2 = n - n.1
  
  CCC = h.wanzheng/n^(-(1/(4+px)))
  h.train = CCC * n.1^(-(1/(4+px)))
  h.test = CCC * n.2^(-(1/(4+px)))
  
  # 30-fold for tunning dynamic thresholding parameter 
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
  
  x = testpoint[shushu,]
  thresholding = seq(threshold.min, threshold.max, by = 0.1)
  
  K = guassionkernel(sweep(sweep(X,2,x,"-"),2,h.wanzheng,"/"))/(prod(h.wanzheng))
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
      XXX = X[a,]
      YYY.res = matrix(1, nrow = n.2, ncol = pp)
      YYY.res = Y[res, ]
      XXX.res = X[res,]
      
      K = guassionkernel(sweep(sweep(XXX,2,x,"-"),2,h.train,"/"))/(prod(h.train))
      
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
      
      
      K.res = guassionkernel(sweep(sweep(XXX.res,2,x,"-"),2,h.test,"/"))/(prod(h.test))
      
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
  p1 <- 5 #10
  # Dimension of Y
  p2 <- 100 #200,300 
  
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
  set.seed(1) # seed1 for models 1,2; seed28 for models 3,4
  testpoint<- matrix(runif(ns * p1)*2-1, ns,p1)
  ##########################################################################################################
  ##########################################################################################################
  
  thresh1 = numeric(ns)
  
  results = bplapply(seq(1:ns), dcm, testpoint, X, Y ,BPPARAM = mcparam)
  
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
  
  # Only used for models 3 and 4
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




