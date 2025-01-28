# Kernel function
guassionkernel = function(x) exp(-0.5 * (x^2))/ sqrt(2 * pi)

# Main funtion
dcm = function(kks1,YY, XX, threshold.min, threshold.max) {
  guassionkernel = function(x) exp(-0.5 * (x^2))/ sqrt(2 * pi)
  library(MASS)
  
  X = XX[kks1:(kks1+99)]
  Y = YY[kks1:(kks1+99),]
  Xs = XX[kks1+100]
  test_Y = YY[kks1+100,]
  
  
  n  = length(X)
  pp = length(Y[1, ])
  error.NW = Y
  sigma.NW.jihe = list()
  
  ######################################
  #### Tuning the bandwidth parameter
  ######################################
  h1 = n^(-0.2) * median(abs(X - median(X)))/0.6745
  bandwidth = seq(h1, h1 + 0.9, by = 0.1)
  
  NN = 30
  pppp = round(pp/200)
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
      nnn1=n
      a = bandwidthA[j, ]  # sample(1:pp, pppp, replace = FALSE, prob = NULL)
      
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
  ### Calculate the estimated covariance matrices for each test point, we only use the soft threshold here.
  ##################################################################################################################
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
  

  x = Xs
  thresholding = seq(threshold.min, threshold.max, by = 0.01)
  error.thresh.soft = numeric(length(thresholding))
    
  for (thresh in 1:length(thresholding)) {
      
    s = thresholding[thresh]
    sum.NWs.soft = 0
      
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
      sigma.NW.jihe.soft1 = matrix(0, ncol= pp, nrow = pp)
      sigma.NW.jihe.res1 = matrix(0, ncol = pp, nrow = pp)
        
      for (i in 1:n.1) {
        sigma.NW.jihe1 = sigma.NW.jihe1 + K[i] * (t(t(YYY[i, ])) %*% t(YYY[i, ]))
      }
      sigma.NW.jihe1 = sigma.NW.jihe1/sum(K)
      sigma.NW.jihe.soft1 = sign(sigma.NW.jihe1) *
        (abs(sigma.NW.jihe1) - s) *
        ((abs(sigma.NW.jihe1) - s) > 0)
      K.res = guassionkernel((XXX.res - x)/h.test)/h.test
        
      for (i in 1:n.2) {
        sigma.NW.jihe.res1 = sigma.NW.jihe.res1 + K.res[i] * (t(t(YYY.res[i,])) %*% t(YYY.res[i, ]))
      }
        
      sigma.NW.jihe.res1 = sigma.NW.jihe.res1/sum(K.res)
      HHH = sigma.NW.jihe.soft1 - sigma.NW.jihe.res1
      sum.NWs.soft = sum.NWs.soft + sum(diag(HHH %*% HHH))
        
    }
      
    error.thresh.soft[thresh] = sum.NWs.soft/n
      
  }
    
  s.soft = thresholding[which.min(error.thresh.soft)]
  sigma.NW.jihe.soft = sign(sigma.NW.jihe) * (abs(sigma.NW.jihe) - s.soft) * ((abs(sigma.NW.jihe) - s.soft) > 0)
  diag(sigma.NW.jihe.soft)=diag_x
  sigma.NW.jihe.soft_C=sigma.NW.jihe.soft+((min(eigen(sigma.NW.jihe.soft)$val)<=0)*(-min(eigen(sigma.NW.jihe.soft)$val)+0.1))*diag(pp)
  
  # Result for Kernel
  w_t = (1/(t(rep(1, pp))%*%solve(sigma.NW.jihe.soft)%*%rep(1, pp)))[1,1]*(solve(sigma.NW.jihe.soft)%*%rep(1, pp))
  r_i_DCM_soft = sum(test_Y*w_t)
  
  # Result for Modified-Kernel
  w_t_c = (1/(t(rep(1, pp))%*%solve(sigma.NW.jihe.soft_C)%*%rep(1, pp)))[1,1]*(solve(sigma.NW.jihe.soft_C)%*%rep(1, pp))
  r_i_DCM_soft_c = sum(test_Y*w_t_c)
  return(list(r_i_DCM_soft,r_i_DCM_soft_c))
  
}

library(readxl)
library(BiocParallel)

# Number of cores
ncores = 40
mcparam = SnowParam(workers = ncores)

# Read data
data <- read_excel("real_data_2.xlsx",col_names = FALSE)

matrix_data <- as.matrix(data)

Y = matrix_data[, 1:(ncol(matrix_data) - 305)]

# Mkt-RF as covariate
X = matrix_data[, (ncol(matrix_data) - 4):(ncol(matrix_data) - 4)]

B=bplapply(seq(1:2416),dcm,Y,X,0.1,0.6,BPPARAM = mcparam)

# Calculate AVR,STD and IR
ans<-matrix(unlist(B), ncol= 2 , byrow= TRUE )
col_means <- apply(ans, 2, mean)
col_sd <- apply(ans, 2, sd)

AVR <- col_means*252
STD <- col_sd*sqrt(252)
IR <- AVR/STD

# Print the result
AVR
STD
IR


