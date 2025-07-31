library(readxl)
library(grf)
library(parallel)
library(foreach)
library(iterators)
library(doParallel)
library(BiocParallel)

# Read data
data <- read_excel("real_data_2.xlsx",col_names = FALSE)

# Sample size
n = 100
n.1 = round(0.7*n)
n.2 = n - n.1
# Dimension of Y
p=200

matrix_data <- as.matrix(data)

Y = matrix_data[, 1:(ncol(matrix_data) - 305)]

X = matrix_data[, (ncol(matrix_data) - 4):ncol(matrix_data)]

out_array_B_new=list()

# Number of cores
ncores = 40
mcparam = SnowParam(workers = ncores)

train_Y = Y[1:100,]
train_X = X[1:100,]
train_Y2<- matrix(0,n,(p*p+p)/2)

for(i in 1:100){
  B<-as.vector(upper.tri(train_Y[i,]%*%t(train_Y[i,]), diag=TRUE))
  C<-as.vector(train_Y[i,]%*%t(train_Y[i,]))
  train_Y2[i,]=C[B]
}


R.forest <-multi_regression_forest(train_X,train_Y,num.trees = 2000)
R2.forest <-multi_regression_forest(train_X,train_Y2,num.trees = 2000)

iter_function = function(kkk,X,Y,R.forest,R2.forest,train_Y,train_X,train_Y2)
{
  
  library(grf)
  
  # 30-fold for tunning dynamic thresholding parameter
  n = 100
  n.1 = round(0.7*n)
  n.2 = n - n.1
  p=200
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
  thresholding = seq(0, 2, by = 0.1)
  
  x = X[kkk,]
  
  EX=predict(R.forest,t(as.matrix(x)))$predictions
  EX2=predict(R2.forest,t(as.matrix(x)))$predictions
  EX_2=matrix(0,p,p)
  EX_2[upper.tri(EX_2, diag=TRUE)] <- EX2
  EX_2[lower.tri(EX_2)] <- t(EX_2)[lower.tri(EX_2)]
  sigma.NW.jihe=EX_2-t(EX)%*%EX
  
  error.thresh.soft = numeric(length(thresholding))
  
  for (thresh in 1:length(thresholding)) {
    s = thresholding[thresh]
    sum.NWs.soft = 0
    for (j in 1:N1) {
      a = TAa[j, ]
      res = TAa.res[j, ]
      YYY = train_Y[a,]
      XXX = train_X[a,]
      YYY2 = train_Y2[a,]
      
      YYY.res = train_Y[res, ]
      XXX.res = train_X[res,]
      YYY2.res = train_Y2[res,]
      
      r.forest <-multi_regression_forest(as.matrix(XXX),as.matrix(YYY),num.trees = 200)
      r2.forest <-multi_regression_forest(as.matrix(XXX),as.matrix(YYY2),num.trees = 200)
      
      EX=predict(r.forest,t(as.matrix(x)))$predictions
      EX2=predict(r2.forest,t(as.matrix(x)))$predictions
      
      EX_2=matrix(0,p,p)
      EX_2[upper.tri(EX_2, diag=TRUE)] <- EX2
      EX_2[lower.tri(EX_2)] <- t(EX_2)[lower.tri(EX_2)]
      sigma.NW.jihe1=EX_2-t(EX)%*%EX
      
      diag_x=diag(sigma.NW.jihe1)
      
      sigma.NW.jihe.soft1 = sign(sigma.NW.jihe1) *
        (abs(sigma.NW.jihe1) - s) *
        ((abs(sigma.NW.jihe1) - s) > 0)
      diag(sigma.NW.jihe.soft1)=diag_x
      
      res.forest <-multi_regression_forest(as.matrix(XXX.res),as.matrix(YYY.res),num.trees = 200)
      res2.forest <-multi_regression_forest(as.matrix(XXX.res),as.matrix(YYY2.res),num.trees = 200)
      
      EX=predict(res.forest,t(as.matrix(x)))$predictions
      EX2=predict(res2.forest,t(as.matrix(x)))$predictions
      
      EX_2=matrix(0,p,p)
      EX_2[upper.tri(EX_2, diag=TRUE)] <- EX2
      EX_2[lower.tri(EX_2)] <- t(EX_2)[lower.tri(EX_2)]
      sigma.NW.jihe.res1=EX_2-t(EX)%*%EX
      
      diag(sigma.NW.jihe.res1)=diag_x
      
      
      HHH = sigma.NW.jihe.soft1 - sigma.NW.jihe.res1
      sum.NWs.soft = sum.NWs.soft + sum(diag(HHH %*% HHH))
      
      
    }
    error.thresh.soft[thresh] = sum.NWs.soft/n
    
  }
  rf.soft = thresholding[which.min(error.thresh.soft)]
  
  diag_x=diag(sigma.NW.jihe)
  
  sigma.NW.jihe.soft = sign(sigma.NW.jihe) * (abs(sigma.NW.jihe) - rf.soft) * ((abs(sigma.NW.jihe) - rf.soft) > 0)
  diag(sigma.NW.jihe.soft)=diag_x
  
  
  sigma.NW.jihe.soft_corr=sigma.NW.jihe.soft+((min(eigen(sigma.NW.jihe.soft)$val)<=0)*(-min(eigen(sigma.NW.jihe.soft)$val)+0.1))*diag(p)
  
  # Static
  m=apply(train_Y,2,mean)
  m_2=(t(train_Y)%*%train_Y)/n
  sigma.jihe = m_2-m%*%t(m)
  
  error.thresh.soft = numeric(length(thresholding))
  
  for (thresh in 1:length(thresholding)) {
    s = thresholding[thresh]
    sum.soft = 0
    for (j in 1:N1) {
      a = TAa[j, ]
      res = TAa.res[j, ]
      YYY = train_Y[a,]
      XXX = train_X[a,]
      
      YYY.res = train_Y[res, ]
      XXX.res = train_X[res,]
      
      m=apply(YYY,2,mean)
      m_2=(t(YYY)%*%YYY)/n.1
      sigma.jihe1 =m_2-m%*%t(m)
      diag_x=diag(sigma.jihe1)
      
      sigma.jihe.soft1 = sign(sigma.jihe1) *(abs(sigma.jihe1) - s) *((abs(sigma.jihe1) - s) > 0)
      diag(sigma.jihe.soft1)=diag_x
      
      
      
      m=apply(YYY.res,2,mean)
      m_2=(t(YYY.res)%*%YYY.res)/n.2
      sigma.jihe.res1 =m_2-m%*%t(m)
      
      HHH = sigma.jihe.soft1 - sigma.jihe.res1
      sum.soft = sum.soft + sum(diag(HHH %*% HHH))
      
    }
    error.thresh.soft[thresh] = sum.soft/n
  }
  diag_x=diag(sigma.jihe)
  
  s.soft = thresholding[which.min(error.thresh.soft)]
  
  
  sigma.jihe.soft = sign(sigma.jihe) * (abs(sigma.jihe) - s.soft) * ((abs(sigma.jihe) - s.soft) > 0)
  diag(sigma.jihe.soft)=diag_x
  
  sigma.jihe.soft_corr=sigma.jihe.soft+(-(min(eigen(sigma.jihe.soft)$val)<0)*min(eigen(sigma.jihe.soft)$val)+0.1)*diag(p)
  
  return(list(sigma.NW.jihe,sigma.NW.jihe.soft_corr,sigma.NW.jihe.soft,sigma.jihe.soft_corr))
}
out_array_B_new <-bplapply(seq(1:100), iter_function,X,Y,R.forest,R2.forest,train_Y,train_X,train_Y2,BPPARAM = mcparam)

