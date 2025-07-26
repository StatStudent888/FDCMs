library(grf)
library(MASS)
library(lava)
library(BiocParallel)
library(parallel)



# Number of cores
ncores = 25 
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

# Sample size
n <- 100 
# Dimension of U
p1 <- 5 #10
# Dimension of Y
p2 <- 100 #200,300

iter_function=function(shu,n,p1,p2){
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
  
  # Generate data
  set.seed(shu)
  X <- matrix(runif(n * p1)*2-1, n, p1)
  Y<- matrix(0,n,p2)
  for (i in 1:n){
    set.seed(i)
    Y[i,]<-mvrnorm(1, Dynamic_mean(X[i,],p2), Dynamic_Covariance2(X[i,],p2))
  }
  Y2<- matrix(0,n,(p2*p2+p2)/2)
  for(i in 1:n){
    B<-as.vector(upper.tri(Y[i,]%*%t(Y[i,]), diag=TRUE))
    C<-as.vector(Y[i,]%*%t(Y[i,]))
    Y2[i,]=C[B]
  }
  
  # train random forest
  R.forest <-multi_regression_forest(X,Y,num.trees = 500) #2000
  R2.forest <-multi_regression_forest(X,Y2,num.trees = 500) #2000
  
  # tunning threshold
  threshold.min=0
  thresholding = seq(threshold.min, 2, by = 0.1)
  threshlength <- length(thresholding)
  # The parameter used to do positive definiteness correction (Modified-Kernel)
  c0=0.001
  
  # Generate test points
  ns=30
  set.seed(1) #seed1 for models 1,2; seed28 for models 3,4
  testpoint<- matrix(runif(ns * p1)*2-1, ns,p1)
  ##########################################################################################################
  ##########################################################################################################
  n.1 = round(0.7*n)
  n.2 = n - n.1
  
  # 30-fold for tunning dynamic thresholding parameter 
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
  
  error.thresh.soft = matrix(0, ncol = ns, nrow = threshlength)
  error.thresh.hard = matrix(0, ncol = ns, nrow = threshlength)
  error.thresh.scad = matrix(0, ncol = ns, nrow = threshlength)
  error.thresh.lasso = matrix(0, ncol = ns, nrow = threshlength)
  
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
    
    res.forest <-multi_regression_forest(as.matrix(XXX.res),as.matrix(YYY.res),num.trees = 200)
    res2.forest <-multi_regression_forest(as.matrix(XXX.res),as.matrix(YYY2.res),num.trees = 200)
    
    for (thresh in 1:threshlength){
      s = thresholding[thresh]
      for (j in 1:ns){
        x = testpoint[j,]
        
        EX=predict(R.forest,t(as.matrix(x)))$predictions
        EX2=predict(R2.forest,t(as.matrix(x)))$predictions
        EX_2=matrix(0,p2,p2)
        EX_2[upper.tri(EX_2, diag=TRUE)] <- EX2
        EX_2[lower.tri(EX_2)] <- t(EX_2)[lower.tri(EX_2)]
        sigma.NW.jihe=EX_2-t(EX)%*%EX
        
        diag_x=diag(sigma.NW.jihe)
        
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
        
        EX=predict(res.forest,t(as.matrix(x)))$predictions
        EX2=predict(res2.forest,t(as.matrix(x)))$predictions
        EX_2=matrix(0,p2,p2)
        EX_2[upper.tri(EX_2, diag=TRUE)] <- EX2
        EX_2[lower.tri(EX_2)] <- t(EX_2)[lower.tri(EX_2)]
        sigma.NW.jihe.res1=EX_2-t(EX)%*%EX
        
        
        diag(sigma.NW.jihe.res1)=diag_x
        
        HHH = sigma.NW.jihe.soft1 - sigma.NW.jihe.res1
        error.thresh.soft[thresh,j] <- error.thresh.soft[thresh,j]+(sum(diag(HHH %*% HHH))/N1)
        
        HHH = sigma.NW.jihe.hard1 - sigma.NW.jihe.res1
        error.thresh.hard[thresh,j] <- error.thresh.hard[thresh,j]+(sum(diag(HHH %*% HHH))/N1)
        
        HHH = sigma.NW.jihe.scad1 - sigma.NW.jihe.res1
        error.thresh.scad[thresh,j] <- error.thresh.scad[thresh,j]+(sum(diag(HHH %*% HHH))/N1)
        
        HHH = sigma.NW.jihe.lasso1 - sigma.NW.jihe.res1
        error.thresh.lasso[thresh,j] <- error.thresh.lasso[thresh,j]+(sum(diag(HHH %*% HHH))/N1)
      }
    }
    rm(r.forest,r2.forest,res.forest,res2.forest)
    gc()
  }
  
  s.soft.index <- apply(error.thresh.soft,2,which.min)
  s.hard.index <- apply(error.thresh.hard,2,which.min)
  s.scad.index <- apply(error.thresh.scad,2,which.min)
  s.lasso.index <- apply(error.thresh.lasso,2,which.min)
  
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
  
  for(kk in 1:ns){
    x = testpoint[kk,]
    
    EX=predict(R.forest,t(as.matrix(x)))$predictions
    EX2=predict(R2.forest,t(as.matrix(x)))$predictions
    EX_2=matrix(0,p2,p2)
    EX_2[upper.tri(EX_2, diag=TRUE)] <- EX2
    EX_2[lower.tri(EX_2)] <- t(EX_2)[lower.tri(EX_2)]
    sigma.NW.jihe=EX_2-t(EX)%*%EX
    diag_x=diag(sigma.NW.jihe)
    
    s.soft = thresholding[s.soft.index[kk]]
    s.hard = thresholding[s.hard.index[kk]]
    s.scad = thresholding[s.scad.index[kk]]
    s.lasso = thresholding[s.lasso.index[kk]]
    
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
    
    
    truth <- Dynamic_Covariance2(testpoint[kk,],p2)
    eig <- eigen(sigma.NW.jihe.soft-truth)
    eig_vals <- eig$values
    spectral_radius_soft[1,kk] <- max(abs(eig_vals))
    F_norm_soft[1,kk]<-sqrt(tr((sigma.NW.jihe.soft-truth)
                               %*%(sigma.NW.jihe.soft-truth)))
    
    
    eig <- eigen(sigma.NW.jihe.hard-truth)
    eig_vals <- eig$values
    spectral_radius_hard[1,kk] <- max(abs(eig_vals))
    F_norm_hard[1,kk]<-sqrt(tr((sigma.NW.jihe.hard-truth)
                               %*%(sigma.NW.jihe.hard-truth)))
    
    eig <- eigen(sigma.NW.jihe.scad-truth)
    eig_vals <- eig$values
    spectral_radius_scad[1,kk] <- max(abs(eig_vals))
    F_norm_scad[1,kk]<-sqrt(tr((sigma.NW.jihe.scad-truth)
                               %*%(sigma.NW.jihe.scad-truth)))
    
    eig <- eigen(sigma.NW.jihe.lasso-truth)
    eig_vals <- eig$values
    spectral_radius_lasso[1,kk] <- max(abs(eig_vals))
    F_norm_lasso[1,kk]<-sqrt(tr((sigma.NW.jihe.lasso-truth)
                                %*%(sigma.NW.jihe.lasso-truth)))
    
    
    
    eig <- eigen(sigma.NW.jihe.soft_c-truth)
    eig_vals <- eig$values
    spectral_radius_soft_1[1,kk] <- max(abs(eig_vals))
    F_norm_soft_1[1,kk]<-sqrt(tr((sigma.NW.jihe.soft_c-truth)
                                 %*%(sigma.NW.jihe.soft_c-truth)))
    
    
    eig <- eigen(sigma.NW.jihe.hard_c-truth)
    eig_vals <- eig$values
    spectral_radius_hard_1[1,kk] <- max(abs(eig_vals))
    F_norm_hard_1[1,kk]<-sqrt(tr((sigma.NW.jihe.hard_c-truth)
                                 %*%(sigma.NW.jihe.hard_c-truth)))
    
    eig <- eigen(sigma.NW.jihe.scad_c-truth)
    eig_vals <- eig$values
    spectral_radius_scad_1[1,kk] <- max(abs(eig_vals))
    F_norm_scad_1[1,kk]<-sqrt(tr((sigma.NW.jihe.scad_c-truth)
                                 %*%(sigma.NW.jihe.scad_c-truth)))
    
    eig <- eigen(sigma.NW.jihe.lasso_c-truth)
    eig_vals <- eig$values
    spectral_radius_lasso_1[1,kk] <- max(abs(eig_vals))
    F_norm_lasso_1[1,kk]<-sqrt(tr((sigma.NW.jihe.lasso_c-truth)
                                  %*%(sigma.NW.jihe.lasso_c-truth)))
    
    FPR_loss_soft[1,kk]<-sum(sigma.NW.jihe.soft!=0 & truth ==0)/sum(truth ==0)
    FPR_loss_hard[1,kk]<-sum(sigma.NW.jihe.hard!=0 & truth ==0)/sum(truth ==0)
    FPR_loss_scad[1,kk]<-sum(sigma.NW.jihe.scad!=0 & truth ==0)/sum(truth ==0)
    FPR_loss_lasso[1,kk]<-sum(sigma.NW.jihe.lasso!=0 & truth ==0)/sum(truth ==0)
    
    TPR_loss_soft[1,kk]<-sum(sigma.NW.jihe.soft!=0 & truth !=0)/sum(truth !=0)
    TPR_loss_hard[1,kk]<-sum(sigma.NW.jihe.hard!=0 & truth !=0)/sum(truth !=0)
    TPR_loss_scad[1,kk]<-sum(sigma.NW.jihe.scad!=0 & truth !=0)/sum(truth !=0)
    TPR_loss_lasso[1,kk]<-sum(sigma.NW.jihe.lasso!=0 & truth !=0)/sum(truth !=0)
  }
  MFL_softk=median(F_norm_soft)
  MSL_softk=median(spectral_radius_soft)
  
  MFL_hardk=median(F_norm_hard)
  MSL_hardk=median(spectral_radius_hard)
  
  MFL_scadk=median(F_norm_scad)
  MSL_scadk=median(spectral_radius_scad)
  
  MFL_lassok=median(F_norm_lasso)
  MSL_lassok=median(spectral_radius_lasso)
  
  
  MFL_soft_ck=median(F_norm_soft_1)
  MSL_soft_ck=median(spectral_radius_soft_1)
  
  MFL_hard_ck=median(F_norm_hard_1)
  MSL_hard_ck=median(spectral_radius_hard_1)
  
  MFL_scad_ck=median(F_norm_scad_1)
  MSL_scad_ck=median(spectral_radius_scad_1)
  
  MFL_lasso_ck=median(F_norm_lasso_1)
  MSL_lasso_ck=median(spectral_radius_lasso_1)
  
  TPR_softk=median(TPR_loss_soft)
  FPR_softk=median(FPR_loss_soft)
  
  TPR_hardk=median(TPR_loss_hard)
  FPR_hardk=median(FPR_loss_hard)
  
  TPR_scadk=median(TPR_loss_scad)
  FPR_scadk=median(FPR_loss_scad)
  
  TPR_lassok=median(TPR_loss_lasso)
  FPR_lassok=median(FPR_loss_lasso)
  
  rm(R.forest,R2.forest)
  gc()
  
  ################################################################################################################################################
  ################################################################################################################################################
  # Static
  # 30-fold for tunning dynamic thresholding parameter 
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
    error.thresh.soft[thresh] = sum.soft/N1
    error.thresh.hard[thresh] = sum.hard/N1
    error.thresh.scad[thresh] = sum.scad/N1
    error.thresh.lasso[thresh] = sum.lasso/N1
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
  
  
  for(kk in 1:ns){ 
    truth <- Dynamic_Covariance2(testpoint[kk,],p2)
    eig <- eigen(sigma.jihe.soft-truth)
    eig_vals <- eig$values
    spectral_radius_s_soft[1,kk] <- max(abs(eig_vals))
    F_norm_s_soft[1,kk]<-sqrt(tr((sigma.jihe.soft-truth)
                                 %*%(sigma.jihe.soft-truth)))
    
    
    eig <- eigen(sigma.jihe.hard-truth)
    eig_vals <- eig$values
    spectral_radius_s_hard[1,kk] <- max(abs(eig_vals))
    F_norm_s_hard[1,kk]<-sqrt(tr(((sigma.jihe.hard-truth)
                                  %*%(sigma.jihe.hard-truth))))
    
    eig <- eigen(sigma.jihe.scad-truth)
    eig_vals <- eig$values
    spectral_radius_s_scad[1,kk] <- max(abs(eig_vals))
    F_norm_s_scad[1,kk]<-sqrt(tr((sigma.jihe.scad-truth)
                                 %*%(sigma.jihe.scad-truth)))
    
    
    eig <- eigen(sigma.jihe.lasso-truth)
    eig_vals <- eig$values
    spectral_radius_s_lasso[1,kk] <- max(abs(eig_vals))
    F_norm_s_lasso[1,kk]<-sqrt(tr(((sigma.jihe.lasso-truth)
                                   %*%(sigma.jihe.lasso-truth))))
    
    FPR_loss_s_soft[1,kk]<-sum(sigma.jihe.soft!=0 & truth ==0)/sum(truth ==0)
    FPR_loss_s_hard[1,kk]<-sum(sigma.jihe.hard!=0 & truth ==0)/sum(truth ==0)
    FPR_loss_s_scad[1,kk]<-sum(sigma.jihe.scad!=0 & truth ==0)/sum(truth ==0)
    FPR_loss_s_lasso[1,kk]<-sum(sigma.jihe.lasso!=0 & truth ==0)/sum(truth ==0)
    
    TPR_loss_s_soft[1,kk]<-sum(sigma.jihe.soft!=0 & truth !=0)/sum(truth !=0)
    TPR_loss_s_hard[1,kk]<-sum(sigma.jihe.hard!=0 & truth !=0)/sum(truth !=0)
    TPR_loss_s_scad[1,kk]<-sum(sigma.jihe.scad!=0 & truth !=0)/sum(truth !=0)
    TPR_loss_s_lasso[1,kk]<-sum(sigma.jihe.lasso!=0 & truth !=0)/sum(truth !=0)
  }
  
  MFL_soft_sk=median(F_norm_s_soft)
  MSL_soft_sk=median(spectral_radius_s_soft)
  
  MFL_hard_sk=median(F_norm_s_hard)
  MSL_hard_sk=median(spectral_radius_s_hard)
  
  
  MFL_scad_sk=median(F_norm_s_scad)
  MSL_scad_sk=median(spectral_radius_s_scad)
  
  MFL_lasso_sk=median(F_norm_s_lasso)
  MSL_lasso_sk=median(spectral_radius_s_lasso)
  
  TPR_soft_sk=median(TPR_loss_s_soft)
  FPR_soft_sk=median(FPR_loss_s_soft)
  
  TPR_hard_sk=median(TPR_loss_s_hard)
  FPR_hard_sk=median(FPR_loss_s_hard)
  
  
  TPR_scad_sk=median(TPR_loss_s_scad)
  FPR_scad_sk=median(FPR_loss_s_scad)
  
  TPR_lasso_sk=median(TPR_loss_s_lasso)
  FPR_lasso_sk=median(FPR_loss_s_lasso)
  
  return(list(MFL_softk,MSL_softk,MFL_hardk,MSL_hardk,MFL_scadk,MSL_scadk,MFL_lassok,MSL_lassok,MFL_soft_ck,MSL_soft_ck,MFL_hard_ck,MSL_hard_ck,
              MFL_scad_ck,MSL_scad_ck,MFL_lasso_ck,MSL_lasso_ck,TPR_softk,FPR_softk,TPR_hardk,FPR_hardk,TPR_scadk,FPR_scadk,TPR_lassok,FPR_lassok,
              MFL_soft_sk,MSL_soft_sk,MFL_hard_sk,MSL_hard_sk,MFL_scad_sk,MSL_scad_sk,MFL_lasso_sk,MSL_lasso_sk,TPR_soft_sk,FPR_soft_sk,TPR_hard_sk,
              FPR_hard_sk,TPR_scad_sk,FPR_scad_sk,TPR_lasso_sk,FPR_lasso_sk))
}

results <- bplapply(seq(1:NN), iter_function,n,p1,p2,BPPARAM = mcparam)

for(ii in 1:NN){
  MFL_soft[ii]=results[[ii]][[1]]
  MSL_soft[ii]=results[[ii]][[2]]
  
  MFL_hard[ii]=results[[ii]][[3]]
  MSL_hard[ii]=results[[ii]][[4]]
  
  MFL_scad[ii]=results[[ii]][[5]]
  MSL_scad[ii]=results[[ii]][[6]]
  
  MFL_lasso[ii]=results[[ii]][[7]]
  MSL_lasso[ii]=results[[ii]][[8]]
  
  
  MFL_soft_c[ii]=results[[ii]][[9]]
  MSL_soft_c[ii]=results[[ii]][[10]]
  
  MFL_hard_c[ii]=results[[ii]][[11]]
  MSL_hard_c[ii]=results[[ii]][[12]]
  
  MFL_scad_c[ii]=results[[ii]][[13]]
  MSL_scad_c[ii]=results[[ii]][[14]]
  
  MFL_lasso_c[ii]=results[[ii]][[15]]
  MSL_lasso_c[ii]=results[[ii]][[16]]
  
  TPR_soft[ii]=results[[ii]][[17]]
  FPR_soft[ii]=results[[ii]][[18]]
  
  TPR_hard[ii]=results[[ii]][[19]]
  FPR_hard[ii]=results[[ii]][[20]]
  
  TPR_scad[ii]=results[[ii]][[21]]
  FPR_scad[ii]=results[[ii]][[22]]
  
  TPR_lasso[ii]=results[[ii]][[23]]
  FPR_lasso[ii]=results[[ii]][[24]]
  
  MFL_soft_s[ii]=results[[ii]][[25]]
  MSL_soft_s[ii]=results[[ii]][[26]]
  
  MFL_hard_s[ii]=results[[ii]][[27]]
  MSL_hard_s[ii]=results[[ii]][[28]]
  
  
  MFL_scad_s[ii]=results[[ii]][[29]]
  MSL_scad_s[ii]=results[[ii]][[30]]
  
  MFL_lasso_s[ii]=results[[ii]][[31]]
  MSL_lasso_s[ii]=results[[ii]][[32]]
  
  TPR_soft_s[ii]=results[[ii]][[33]]
  FPR_soft_s[ii]=results[[ii]][[34]]
  
  TPR_hard_s[ii]=results[[ii]][[35]]
  FPR_hard_s[ii]=results[[ii]][[36]]
  
  
  TPR_scad_s[ii]=results[[ii]][[37]]
  FPR_scad_s[ii]=results[[ii]][[38]]
  
  TPR_lasso_s[ii]=results[[ii]][[39]]
  FPR_lasso_s[ii]=results[[ii]][[40]]
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

