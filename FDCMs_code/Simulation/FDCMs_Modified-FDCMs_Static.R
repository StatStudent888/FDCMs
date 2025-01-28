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

library(grf)
library(MASS)
library(lava)
library(BiocParallel)
library(parallel)

# Number of cores
ncores = 40 
mcparam = SnowParam(workers = ncores)
# Number of repetitions
NN=50 

# MFL,MSL for FDCMs
MFL_soft=matrix(0, ncol = NN, nrow = 1)
MSL_soft=matrix(0, ncol = NN, nrow = 1)

MFL_hard=matrix(0, ncol = NN, nrow = 1)
MSL_hard=matrix(0, ncol = NN, nrow = 1)

MFL_lasso=matrix(0, ncol = NN, nrow = 1)
MSL_lasso=matrix(0, ncol = NN, nrow = 1)

MFL_scad=matrix(0, ncol = NN, nrow = 1)
MSL_scad=matrix(0, ncol = NN, nrow = 1)


# MFL,MSL for Modified-FDCMs 
MFL_soft_c=matrix(0, ncol = NN, nrow = 1)
MSL_soft_c=matrix(0, ncol = NN, nrow = 1)

MFL_hard_c=matrix(0, ncol = NN, nrow = 1)
MSL_hard_c=matrix(0, ncol = NN, nrow = 1)

MFL_lasso_c=matrix(0, ncol = NN, nrow = 1)
MSL_lasso_c=matrix(0, ncol = NN, nrow = 1)

MFL_scad_c=matrix(0, ncol = NN, nrow = 1)
MSL_scad_c=matrix(0, ncol = NN, nrow = 1)


# MFL,MSL for Static
MFL_soft_s=matrix(0, ncol = NN, nrow = 1)
MSL_soft_s=matrix(0, ncol = NN, nrow = 1)

MFL_hard_s=matrix(0, ncol = NN, nrow = 1)
MSL_hard_s=matrix(0, ncol = NN, nrow = 1)

MFL_lasso_s=matrix(0, ncol = NN, nrow = 1)
MSL_lasso_s=matrix(0, ncol = NN, nrow = 1)

MFL_scad_s=matrix(0, ncol = NN, nrow = 1)
MSL_scad_s=matrix(0, ncol = NN, nrow = 1)


# TPR,FPR for FDCMs (only used for model 3 and 4)
TPR_soft=matrix(0, ncol = NN, nrow = 1)
FPR_soft=matrix(0, ncol = NN, nrow = 1)

TPR_hard=matrix(0, ncol = NN, nrow = 1)
FPR_hard=matrix(0, ncol = NN, nrow = 1)

TPR_lasso=matrix(0, ncol = NN, nrow = 1)
FPR_lasso=matrix(0, ncol = NN, nrow = 1)

TPR_scad=matrix(0, ncol = NN, nrow = 1)
FPR_scad=matrix(0, ncol = NN, nrow = 1)


# TPR,FPR for Static (only used for model 3 and 4)
TPR_soft_s=matrix(0, ncol = NN, nrow = 1)
FPR_soft_s=matrix(0, ncol = NN, nrow = 1)

TPR_hard_s=matrix(0, ncol = NN, nrow = 1)
FPR_hard_s=matrix(0, ncol = NN, nrow = 1)

TPR_lasso_s=matrix(0, ncol = NN, nrow = 1)
FPR_lasso_s=matrix(0, ncol = NN, nrow = 1)

TPR_scad_s=matrix(0, ncol = NN, nrow = 1)
FPR_scad_s=matrix(0, ncol = NN, nrow = 1)



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
  Y<- matrix(0,n,p2)
  for (i in 1:n)
  {
    set.seed(i)
    Y[i,]<-mvrnorm(1, Dynamic_mean(X[i,],p2), Dynamic_Covariance2(X[i,],p2))
  }
  Y2<- matrix(0,n,(p2*p2+p2)/2)
  for(i in 1:n){
    B<-as.vector(upper.tri(Y[i,]%*%t(Y[i,]), diag=TRUE))
    C<-as.vector(Y[i,]%*%t(Y[i,]))
    Y2[i,]=C[B]
  }
  
  # Generate test points
  ns=30
  set.seed(ii+1) # used for model 1,4; set.seed(ii+6) for model 2,3
  testpoint<- matrix(runif(ns * p1)*2-1, ns,p1)
  ##########################################################################################################
  ##########################################################################################################
  n.1 = round(0.7*n)
  n.2 = n - n.1
  
  # 30-fold for tuning dynamic thresholding parameter 
  N1 = 30
  TAa = matrix(0, nrow = N1, ncol = n.1)
  TAa.res = matrix(0, nrow = N1, ncol = n.2)
  
  for (j in 1:N1) {
    set.seed(j)
    a = sample(1:n, n.1, replace = FALSE, prob = NULL)
    a = sort(a)
    res = (1:n)[is.na(pmatch(1:n, a))]
    TAa[j, ] = a
    TAa.res[j, ] = res
  }
  ##########################################################################################################
  #FDCMs function
  iter_function = function(shu, testpoint, N1, TAa, thresh1, R.forest, threshold.min,
                           TAa.res, X, Y,  n.1, n.2, n, p1, p2, R2.forest){
    # SCAD threshold function
    scad_penalty <- function(x, lambda, a) {
      abs_x <- abs(x)
      ifelse(abs_x <= 2*lambda, sign(x) *(abs_x-lambda)*((abs_x-lambda)>0), 
             ifelse(abs_x <= (a * lambda), ((a - 1) * x -sign(x)*a* lambda)/ (a-2), x))
    }
    
    x = testpoint[shu,]
    
    EX=predict(R.forest,t(as.matrix(x)))$predictions
    EX2=predict(R2.forest,t(as.matrix(x)))$predictions
    EX_2=matrix(0,p2,p2)
    EX_2[upper.tri(EX_2, diag=TRUE)] <- EX2
    EX_2[lower.tri(EX_2)] <- t(EX_2)[lower.tri(EX_2)]
    sigma.NW.jihe=EX_2-t(EX)%*%EX
    
    diag_x=diag(sigma.NW.jihe)
    
    thresholding = seq(threshold.min, 2, by = 0.1)
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
        YYY = Y[a,]
        XXX = X[a,]
        YYY2 = Y2[a,]
        
        YYY.res = Y[res, ]
        XXX.res = X[res,]
        YYY2.res = Y2[res,]
        
        r.forest <-multi_regression_forest(as.matrix(XXX),as.matrix(YYY),num.trees = 200)
        r2.forest <-multi_regression_forest(as.matrix(XXX),as.matrix(YYY2),num.trees = 200)
        
        EX=predict(r.forest,t(as.matrix(x)))$predictions
        EX2=predict(r2.forest,t(as.matrix(x)))$predictions
        
        EX_2=matrix(0,p2,p2)
        EX_2[upper.tri(EX_2, diag=TRUE)] <- EX2
        EX_2[lower.tri(EX_2)] <- t(EX_2)[lower.tri(EX_2)]
        sigma.NW.jihe1=EX_2-t(EX)%*%EX
        
        
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
        
        res.forest <-multi_regression_forest(as.matrix(XXX.res),as.matrix(YYY.res),num.trees = 200)
        res2.forest <-multi_regression_forest(as.matrix(XXX.res),as.matrix(YYY2.res),num.trees = 200)
        
        EX=predict(res.forest,t(as.matrix(x)))$predictions
        EX2=predict(res2.forest,t(as.matrix(x)))$predictions
        
        EX_2=matrix(0,p2,p2)
        EX_2[upper.tri(EX_2, diag=TRUE)] <- EX2
        EX_2[lower.tri(EX_2)] <- t(EX_2)[lower.tri(EX_2)]
        sigma.NW.jihe.res1=EX_2-t(EX)%*%EX
        
        
        diag(sigma.NW.jihe.res1)=diag_x
        
        
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
    sigma.NW.jihe.soft_c=sigma.NW.jihe.soft+diag(rep(-a,p2)+c0)*(a<0)
    
    sigma.NW.jihe.hard = sigma.NW.jihe *((abs(sigma.NW.jihe) - s.hard) > 0)
    diag(sigma.NW.jihe.hard)=diag_x
    a=min(eigen(sigma.NW.jihe.hard)$values)
    sigma.NW.jihe.hard_c=sigma.NW.jihe.hard+diag(rep(-a,p2)+c0)*(a<0)
    
    sigma.NW.jihe.scad = scad_penalty(sigma.NW.jihe,s.scad,3.7)
    diag(sigma.NW.jihe.scad)=diag_x
    a=min(eigen(sigma.NW.jihe.scad)$values)
    sigma.NW.jihe.scad_c=sigma.NW.jihe.scad+diag(rep(-a,p2)+c0)*(a<0)
    
    sigma.NW.jihe.lasso = sign(sigma.NW.jihe) *(abs(sigma.NW.jihe) - s.lasso^2/abs(sigma.NW.jihe))*((abs(sigma.NW.jihe) - s.lasso^2/abs(sigma.NW.jihe)) > 0)
    diag(sigma.NW.jihe.lasso)=diag_x
    a=min(eigen(sigma.NW.jihe.lasso)$values)
    sigma.NW.jihe.lasso_c=sigma.NW.jihe.lasso+diag(rep(-a,p2)+c0)*(a<0)
    
    return(list(sigma.NW.jihe.soft,sigma.NW.jihe.hard,sigma.NW.jihe.scad,sigma.NW.jihe.lasso,sigma.NW.jihe.soft_c,sigma.NW.jihe.hard_c,sigma.NW.jihe.scad_c,sigma.NW.jihe.lasso_c))
  }
  
  threshold.min=0
  
  thresh1 = numeric(ns)
  
  R.forest <-multi_regression_forest(X,Y)
  R2.forest <-multi_regression_forest(X,Y2)
  
  results = bplapply(seq(1:ns), iter_function,
                     testpoint, N1, TAa,thresh1,R.forest,threshold.min,
                     TAa.res, X, Y,  n.1, n.2, n, p1, p2,R2.forest
                     ,BPPARAM = mcparam)
  
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
  
  
  ################################################################################################################################################
  ################################################################################################################################################
  # Static
  # 30-fold for tuning dynamic thresholding parameter 
  N1=30
  TAA = matrix(0, nrow = N1, ncol = n.1)
  TAA.res = matrix(0, nrow = N1, ncol = n.2)
  
  for (j in 1:N1) {
    set.seed(j)
    a = sample(1:n, n.1, replace = FALSE, prob = NULL)
    a = sort(a)
    res = (1:n)[is.na(pmatch(1:n, a))]
    TAA[j, ] = a
    TAA.res[j, ] = res
  }
  m=apply(Y,2,mean)
  m_2=(t(Y)%*%Y)/n
  sigma.jihe = m_2-m%*%t(m)
  
  
  thresholding = seq(0,2, by = 0.1)
  error.thresh.soft = numeric(length(thresholding))
  error.thresh.hard = numeric(length(thresholding))
  error.thresh.scad = numeric(length(thresholding))
  error.thresh.lasso = numeric(length(thresholding))
  
  for (thresh in 1:length(thresholding)) {
    s = thresholding[thresh]
    sum.soft = 0
    sum.hard = 0
    sum.scad = 0
    sum.lasso = 0
    for (j in 1:N1) {
      a = TAA[j, ]
      res = TAA.res[j, ]
      YYY = Y[a,]
      XXX = X[a,]
      
      YYY.res = Y[res, ]
      XXX.res = X[res,]
      
      sigma.jihe1 = matrix(0, ncol = p2, nrow = p2)
      sigma.jihe.soft1 = matrix(0, ncol= p2, nrow = p2)
      sigma.jihe.res1 = matrix(0, ncol = p2, nrow = p2)
      
      
      m=apply(YYY,2,mean)
      m_2=(t(YYY)%*%YYY)/n.1
      sigma.jihe1 =m_2-m%*%t(m)
      diag_x=diag(sigma.jihe1)
      
      sigma.jihe.soft1 = sign(sigma.jihe1) *(abs(sigma.jihe1) - s) *((abs(sigma.jihe1) - s) > 0)
      diag(sigma.jihe.soft1)=diag_x
      
      sigma.jihe.hard1 = (sigma.jihe1) *((abs(sigma.jihe1) - s) > 0)
      diag(sigma.jihe.hard1)=diag_x
      
      sigma.jihe.scad1 = scad_penalty(sigma.jihe1,s,3.7)
      diag(sigma.jihe.scad1)=diag_x
      
      sigma.jihe.lasso1 = sign(sigma.jihe1) *(abs(sigma.jihe1) - s^2/abs(sigma.jihe1))*((abs(sigma.jihe1) - s^2/abs(sigma.jihe1)) > 0)
      diag(sigma.jihe.lasso1)=diag_x
      
      
      m=apply(YYY.res,2,mean)
      m_2=(t(YYY.res)%*%YYY.res)/n.2
      sigma.jihe.res1 =m_2-m%*%t(m)
      diag(sigma.jihe.res1)=diag_x
      
      HHH = sigma.jihe.soft1 - sigma.jihe.res1
      sum.soft = sum.soft + sum(diag(HHH %*% HHH))
      
      HHH = sigma.jihe.hard1 - sigma.jihe.res1
      sum.hard = sum.hard + sum(diag(HHH %*% HHH))
      
      HHH = sigma.jihe.scad1 - sigma.jihe.res1
      sum.scad = sum.scad + sum(diag(HHH %*% HHH))
      
      HHH = sigma.jihe.lasso1 - sigma.jihe.res1
      sum.lasso = sum.lasso + sum(diag(HHH %*% HHH))
      
    }
    error.thresh.soft[thresh] = sum.soft/n
    error.thresh.hard[thresh] = sum.hard/n
    error.thresh.scad[thresh] = sum.scad/n
    error.thresh.lasso[thresh] = sum.lasso/n
  }
  diag_x=diag(sigma.jihe)
  
  s.soft = thresholding[which.min(error.thresh.soft)]
  
  s.hard = thresholding[which.min(error.thresh.hard)]
  
  s.scad = thresholding[which.min(error.thresh.scad)]
  
  s.lasso = thresholding[which.min(error.thresh.lasso)]
  
  
  sigma.jihe.soft = sign(sigma.jihe) * (abs(sigma.jihe) - s.soft) * ((abs(sigma.jihe) - s.soft) > 0)
  diag(sigma.jihe.soft)=diag_x
  
  sigma.jihe.hard = sigma.jihe *((abs(sigma.jihe) - s.hard) > 0)
  diag(sigma.jihe.hard)=diag_x
  
  sigma.jihe.scad = scad_penalty(sigma.jihe,s.scad,3.7)
  diag(sigma.jihe.scad)=diag_x
  
  sigma.jihe.lasso = sign(sigma.jihe) *(abs(sigma.jihe) - s.lasso^2/abs(sigma.jihe))*((abs(sigma.jihe) - s.lasso^2/abs(sigma.jihe)) > 0)
  diag(sigma.jihe.lasso)=diag_x
  
  F_norm_s_soft<-matrix(0, ncol =ns, nrow = 1)
  F_norm_s_hard<-matrix(0, ncol =ns, nrow = 1)
  F_norm_s_scad<-matrix(0, ncol =ns, nrow = 1)
  F_norm_s_lasso<-matrix(0, ncol =ns, nrow = 1)
  
  spectral_radius_s_soft<-matrix(0, ncol =ns, nrow = 1)
  spectral_radius_s_hard<-matrix(0, ncol =ns, nrow = 1)
  spectral_radius_s_scad<-matrix(0, ncol =ns, nrow = 1)
  spectral_radius_s_lasso<-matrix(0, ncol =ns, nrow = 1)
  
  FPR_loss_s_soft<-matrix(0, ncol =ns, nrow = 1)
  FPR_loss_s_hard<-matrix(0, ncol =ns, nrow = 1)
  FPR_loss_s_scad<-matrix(0, ncol =ns, nrow = 1)
  FPR_loss_s_lasso<-matrix(0, ncol =ns, nrow = 1)
  
  TPR_loss_s_soft<-matrix(0, ncol =ns, nrow = 1)
  TPR_loss_s_hard<-matrix(0, ncol =ns, nrow = 1)
  TPR_loss_s_scad<-matrix(0, ncol =ns, nrow = 1)
  TPR_loss_s_lasso<-matrix(0, ncol =ns, nrow = 1)
  
  
  for(kk in 1:ns)
  { 
    eig <- eigen(sigma.jihe.soft-Dynamic_Covariance2(testpoint[kk,],p2))
    eig_vals <- eig$values
    spectral_radius_s_soft[1,kk] <- max(abs(eig_vals))
    F_norm_s_soft[1,kk]<-sqrt(tr((sigma.jihe.soft-Dynamic_Covariance2(testpoint[kk,],p2))
                                 %*%(sigma.jihe.soft-Dynamic_Covariance2(testpoint[kk,],p2))))
    
    
    eig <- eigen(sigma.jihe.hard-Dynamic_Covariance2(testpoint[kk,],p2))
    eig_vals <- eig$values
    spectral_radius_s_hard[1,kk] <- max(abs(eig_vals))
    F_norm_s_hard[1,kk]<-sqrt(tr(((sigma.jihe.hard-Dynamic_Covariance2(testpoint[kk,],p2))
                                  %*%(sigma.jihe.hard-Dynamic_Covariance2(testpoint[kk,],p2)))))
    
    eig <- eigen(sigma.jihe.scad-Dynamic_Covariance2(testpoint[kk,],p2))
    eig_vals <- eig$values
    spectral_radius_s_scad[1,kk] <- max(abs(eig_vals))
    F_norm_s_scad[1,kk]<-sqrt(tr((sigma.jihe.scad-Dynamic_Covariance2(testpoint[kk,],p2))
                                 %*%(sigma.jihe.scad-Dynamic_Covariance2(testpoint[kk,],p2))))
    
    
    eig <- eigen(sigma.jihe.lasso-Dynamic_Covariance2(testpoint[kk,],p2))
    eig_vals <- eig$values
    spectral_radius_s_lasso[1,kk] <- max(abs(eig_vals))
    F_norm_s_lasso[1,kk]<-sqrt(tr(((sigma.jihe.lasso-Dynamic_Covariance2(testpoint[kk,],p2))
                                   %*%(sigma.jihe.lasso-Dynamic_Covariance2(testpoint[kk,],p2)))))
    
    FPR_loss_s_soft[1,kk]<-sum(sigma.jihe.soft!=0 & Dynamic_Covariance2(testpoint[kk,],p2) ==0)/sum(Dynamic_Covariance2(testpoint[kk,],p2) ==0)
    FPR_loss_s_hard[1,kk]<-sum(sigma.jihe.hard!=0 & Dynamic_Covariance2(testpoint[kk,],p2) ==0)/sum(Dynamic_Covariance2(testpoint[kk,],p2) ==0)
    FPR_loss_s_scad[1,kk]<-sum(sigma.jihe.scad!=0 & Dynamic_Covariance2(testpoint[kk,],p2) ==0)/sum(Dynamic_Covariance2(testpoint[kk,],p2) ==0)
    FPR_loss_s_lasso[1,kk]<-sum(sigma.jihe.lasso!=0 & Dynamic_Covariance2(testpoint[kk,],p2) ==0)/sum(Dynamic_Covariance2(testpoint[kk,],p2) ==0)
    
    TPR_loss_s_soft[1,kk]<-sum(sigma.jihe.soft!=0 & Dynamic_Covariance2(testpoint[kk,],p2) !=0)/sum(Dynamic_Covariance2(testpoint[kk,],p2) !=0)
    TPR_loss_s_hard[1,kk]<-sum(sigma.jihe.hard!=0 & Dynamic_Covariance2(testpoint[kk,],p2) !=0)/sum(Dynamic_Covariance2(testpoint[kk,],p2) !=0)
    TPR_loss_s_scad[1,kk]<-sum(sigma.jihe.scad!=0 & Dynamic_Covariance2(testpoint[kk,],p2) !=0)/sum(Dynamic_Covariance2(testpoint[kk,],p2) !=0)
    TPR_loss_s_lasso[1,kk]<-sum(sigma.jihe.lasso!=0 & Dynamic_Covariance2(testpoint[kk,],p2) !=0)/sum(Dynamic_Covariance2(testpoint[kk,],p2) !=0)
  }
  
  MFL_soft_s[ii]=median(F_norm_s_soft)
  MSL_soft_s[ii]=median(spectral_radius_s_soft)
  
  MFL_hard_s[ii]=median(F_norm_s_hard)
  MSL_hard_s[ii]=median(spectral_radius_s_hard)
  
  
  MFL_scad_s[ii]=median(F_norm_s_scad)
  MSL_scad_s[ii]=median(spectral_radius_s_scad)
  
  MFL_lasso_s[ii]=median(F_norm_s_lasso)
  MSL_lasso_s[ii]=median(spectral_radius_s_lasso)
  
  TPR_soft_s[ii]=median(TPR_loss_s_soft)
  FPR_soft_s[ii]=median(FPR_loss_s_soft)
  
  TPR_hard_s[ii]=median(TPR_loss_s_hard)
  FPR_hard_s[ii]=median(FPR_loss_s_hard)
  
  
  TPR_scad_s[ii]=median(TPR_loss_s_scad)
  FPR_scad_s[ii]=median(FPR_loss_s_scad)
  
  TPR_lasso_s[ii]=median(TPR_loss_s_lasso)
  FPR_lasso_s[ii]=median(FPR_loss_s_lasso)
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

MFL_soft_s_mean <- mean(MFL_soft_s)
MSL_soft_s_mean <- mean(MSL_soft_s)
MFL_hard_s_mean <- mean(MFL_hard_s)
MSL_hard_s_mean <- mean(MSL_hard_s)
MFL_lasso_s_mean <- mean(MFL_lasso_s)
MSL_lasso_s_mean <- mean(MSL_lasso_s)
MFL_scad_s_mean <- mean(MFL_scad_s)
MSL_scad_s_mean <- mean(MSL_scad_s)

TPR_soft_mean <- mean(TPR_soft)
FPR_soft_mean <- mean(FPR_soft)
TPR_hard_mean <- mean(TPR_hard)
FPR_hard_mean <- mean(FPR_hard)
TPR_lasso_mean <- mean(TPR_lasso)
FPR_lasso_mean <- mean(FPR_lasso)
TPR_scad_mean <- mean(TPR_scad)
FPR_scad_mean <- mean(FPR_scad)

TPR_soft_s_mean <- mean(TPR_soft_s)
FPR_soft_s_mean <- mean(FPR_soft_s)
TPR_hard_s_mean <- mean(TPR_hard_s)
FPR_hard_s_mean <- mean(FPR_hard_s)
TPR_lasso_s_mean <- mean(TPR_lasso_s)
FPR_lasso_s_mean <- mean(FPR_lasso_s)
TPR_scad_s_mean <- mean(TPR_scad_s)
FPR_scad_s_mean <- mean(FPR_scad_s)

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

MFL_soft_s_sd <- sd(MFL_soft_s)
MSL_soft_s_sd <- sd(MSL_soft_s)
MFL_hard_s_sd <- sd(MFL_hard_s)
MSL_hard_s_sd <- sd(MSL_hard_s)
MFL_lasso_s_sd <- sd(MFL_lasso_s)
MSL_lasso_s_sd <- sd(MSL_lasso_s)
MFL_scad_s_sd <- sd(MFL_scad_s)
MSL_scad_s_sd <- sd(MSL_scad_s)

TPR_soft_sd <- sd(TPR_soft)
FPR_soft_sd <- sd(FPR_soft)
TPR_hard_sd <- sd(TPR_hard)
FPR_hard_sd <- sd(FPR_hard)
TPR_lasso_sd <- sd(TPR_lasso)
FPR_lasso_sd <- sd(FPR_lasso)
TPR_scad_sd <- sd(TPR_scad)
FPR_scad_sd <- sd(FPR_scad)

TPR_soft_s_sd <- sd(TPR_soft_s)
FPR_soft_s_sd <- sd(FPR_soft_s)
TPR_hard_s_sd <- sd(TPR_hard_s)
FPR_hard_s_sd <- sd(FPR_hard_s)
TPR_lasso_s_sd <- sd(TPR_lasso_s)
FPR_lasso_s_sd <- sd(FPR_lasso_s)
TPR_scad_s_sd <- sd(TPR_scad_s)
FPR_scad_s_sd <- sd(FPR_scad_s)

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

x5 <- paste(sprintf("%0.2f",MFL_hard_s_mean),"(",sprintf("%0.2f",MFL_hard_s_sd),")","&",
            sprintf("%0.2f",MFL_lasso_s_mean),"(",sprintf("%0.2f",MFL_lasso_s_sd),")","&",
            sprintf("%0.2f",MFL_scad_s_mean),"(",sprintf("%0.2f",MFL_scad_s_sd),")","&&",
            sprintf("%0.2f",MFL_soft_s_mean),"(",sprintf("%0.2f",MFL_soft_s_sd),")",sep="")

x6 <- paste(sprintf("%0.2f",MSL_hard_s_mean),"(",sprintf("%0.2f",MSL_hard_s_sd),")","&",
            sprintf("%0.2f",MSL_lasso_s_mean),"(",sprintf("%0.2f",MSL_lasso_s_sd),")","&",
            sprintf("%0.2f",MSL_scad_s_mean),"(",sprintf("%0.2f",MSL_scad_s_sd),")","&&",
            sprintf("%0.2f",MSL_soft_s_mean),"(",sprintf("%0.2f",MSL_soft_s_sd),")",sep="")

x7 <- paste(sprintf("%0.2f",TPR_hard_mean),"(",sprintf("%0.2f",TPR_hard_sd),")","&",
            sprintf("%0.2f",TPR_lasso_mean),"(",sprintf("%0.2f",TPR_lasso_sd),")","&",
            sprintf("%0.2f",TPR_scad_mean),"(",sprintf("%0.2f",TPR_scad_sd),")","&&",
            sprintf("%0.2f",TPR_soft_mean),"(",sprintf("%0.2f",TPR_soft_sd),")",sep="")

x8 <- paste(sprintf("%0.2f",FPR_hard_mean),"(",sprintf("%0.2f",FPR_hard_sd),")","&",
            sprintf("%0.2f",FPR_lasso_mean),"(",sprintf("%0.2f",FPR_lasso_sd),")","&",
            sprintf("%0.2f",FPR_scad_mean),"(",sprintf("%0.2f",FPR_scad_sd),")","&&",
            sprintf("%0.2f",FPR_soft_mean),"(",sprintf("%0.2f",FPR_soft_sd),")",sep="")

x9 <- paste(sprintf("%0.2f",TPR_hard_s_mean),"(",sprintf("%0.2f",TPR_hard_s_sd),")","&",
            sprintf("%0.2f",TPR_lasso_s_mean),"(",sprintf("%0.2f",TPR_lasso_s_sd),")","&",
            sprintf("%0.2f",TPR_scad_s_mean),"(",sprintf("%0.2f",TPR_scad_s_sd),")","&&",
            sprintf("%0.2f",TPR_soft_s_mean),"(",sprintf("%0.2f",TPR_soft_s_sd),")",sep="")

x10 <- paste(sprintf("%0.2f",FPR_hard_s_mean),"(",sprintf("%0.2f",FPR_hard_s_sd),")","&",
            sprintf("%0.2f",FPR_lasso_s_mean),"(",sprintf("%0.2f",FPR_lasso_s_sd),")","&",
            sprintf("%0.2f",FPR_scad_s_mean),"(",sprintf("%0.2f",FPR_scad_s_sd),")","&&",
            sprintf("%0.2f",FPR_soft_s_mean),"(",sprintf("%0.2f",FPR_soft_s_sd),")",sep="")

# Print the result
x1
x2
x3
x4
x5
x6
x7
x8
x9
x10

