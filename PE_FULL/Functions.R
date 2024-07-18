library(glmnet)
#library(ncvreg)
library(MASS)
#library(DescTools)
print("True load")
##############Functions
gen_error<-function(N,p,rho){
  X = matrix(NA,N,p)
  X[,1] = rnorm(N)
  for(ii in 2:p){
    X[,ii] = rho*X[,(ii-1)] + sqrt(1-rho^2)*rnorm(N)
  }
  return(X)
}
trace = function(S)
{
  s = sum(diag(S))
  return(s)
}
#################User functions
Inference_Npe = function(X,Y,CovMatrix,cv,M,seednum)#For n=200,p<=400 moderate dimensional cases. Otherwise can not allocate enough memory.
{
  X = X
  Y = as.vector(Y)
  n = dim(X)[1]
  p = dim(X)[2]
  Omega = Omegahat = Resids = diag(1,p)#Initial storage
  indexX = 1:p
  X1 = X
  Y1 = Y
  n1 = n
  if(cv == "Once")
  {
    n1 = n/2
    X1 = X[1:n1,]
    Y1 = Y[1:n1]   
  }
  if(cv == "Twice")
  {
    n1 = n/2
    X1 = X[-c(1:n1),]
    Y1 = Y[-c(1:n1)]   
    
  }
  if(cv == "Once")
  {
    X = X[-c(1:n1),]
    Y = Y[-c(1:n1)]  
  }
  if(cv == "Twice")
  {
    X = X[c(1:n1),]
    Y = Y[c(1:n1)]  
  }
  if(cv == "Full")
  {
    n1 = n
    X1 = X
    Y1 = Y
  }
  ################Choice of working matrix M
  if(M == "True")
  {
    Omega = solve(CovMatrix)
  }
  if(M == "I")
  {
    Omega = diag(p)
  }
  if(M == "Hat")
  {
    indexX = 1:p
    for(k in 1:p)
    {
      indexZ = indexX[-k]
      Xk = X1[,k]
      Zk = X1[,-k]
      modelx = cv.glmnet(Zk,Xk,nlambda = 10, nfolds = 5)
      lambda = modelx$lambda.min
      hatgamma = coef(modelx, s = "lambda.min")#as.numeric(modelx$beta)
      res = Xk - Zk%*%hatgamma[-1] - hatgamma[1]
      Resids[k,k] = 1/(sum(res^2)/n1+lambda*sum(abs(hatgamma)))
      Omegahat[k,k] = 1
      Omegahat[k,indexZ] = -hatgamma[-1]
    }
    Omega = Resids%*%Omegahat   
  }
  if(M == "I+D1")
  {
    Omega = diag(p) + diag(abs(rnorm(p,0,1)))
  }
  if(M == "I+D2")
  {
    U = rnorm(p,0,1/sqrt(p))
    Omega = diag(p) + U%*%t(U)#solve(CovMatrix) + diag(abs(rnorm(p,0,1)))
  }
  ################
  Z = solve(t(C)%*%Omega%*%C)%*%t(C)%*%Omega%*%t(X)
  W = t(X)-C%*%Z
  mY = mean(Y)
  mZ = apply(Z,1,mean)
  mW = apply(W,1,mean)
  Z1 = Z
  W1 = W
  for(k in 1:n1)
  {
    if(r!=1)
    {
      Z1[,k] = Z[,k]- mZ
    }
    if(r==1)
    {
      Z1[k] = Z[k]-mZ
    }
    W1[,k] = W1[,k] - mW
  }
  modelx = cv.glmnet(X1,Y1,nfolds=5)
  lambda = modelx$lambda.min
  modelx = glmnet(X1,Y1,lambda = lambda)
  beta1 = as.numeric(modelx$beta)
  res = Y-mY-t(Z1)%*%gamma - t(W1)%*%beta1 
  D = sapply(1:ncol(Z),function(x) Z1[,x] * res[x])
  n = dim(D)[2]
  Dcross = apply(D,2,function(x) crossprod(x,D))
  diag(Dcross) = 0
  Tn = sum(Dcross)/(n*(n-1))
  trace_est = sum(Dcross^2)/(n*(n-1))
  estTB2 = 2*(trace_est)/(n*(n-1))
  Tn = Tn/sqrt(estTB2)
  size_npe =  ifelse(abs(Tn)>=qnorm(0.975),1,0)
  list(size_npe = size_npe,stat_npe = Tn,sigma = trace_est)
}

#########################
#########################
#Power enhancement construction under different scenarios
Infer_decor_1 = function(X,Y,seed,spilit)
{
  X = X
  Y = Y
  n = dim(X)[1]
  p = dim(X)[2]
Omega = diag(1,p)
  if(spilit == "Once")
  {
    n1 = n/2
    X1 = X[1:n1,]
    Y1 = Y[1:n1]   
    X = X[-c(1:n1),]
    Y = Y[-c(1:n1)]  
  }
  if(spilit == "Twice")
  {
    n1 = n/2
    X1 = X[-c(1:n1),]
    Y1 = Y[-c(1:n1)]  
    X = X[c(1:n1),]
    Y = Y[c(1:n1)]  
  }
  Z = solve(t(C)%*%Omega%*%C)%*%t(C)%*%Omega%*%t(X)
  W = t(X)-C%*%Z
  mY = mean(Y)
  mZ = apply(Z,1,mean)
  mW = apply(W,1,mean)
  Z1 = Z
  W1 = W
  for(k in 1:n1)
  {
    if(r!=1)
    {
      Z1[,k] = Z[,k]- mZ
    }
    if(r==1)
    {
      Z1[k] = Z[k]-mZ
    }
    W1[,k] = W1[,k] - mW
  }
  res = Y-mY-t(Z1)%*%gamma
  modelx = cv.glmnet(X1,Y1,nfolds=5)
  lambda = modelx$lambda.min
  modelx = glmnet(X1,Y1,lambda = lambda)
  a0 = modelx$a0
  beta1 = as.numeric(modelx$beta)
  res = Y-a0-t(Z1)%*%gamma - t(W1)%*%beta1 
  Score_stat = Pvalue = NULL
  for(kk in 1:r)
  {
    time0 = Sys.time()
    ##########################
    if(model==1|model==2)
    {
      ck = C[,kk]
      Pk = ck%*%solve(t(ck)%*%ck)%*%t(ck)
      Ak = solve(t(ck)%*%ck)%*%t(ck)
      Zk = Ak%*%t(X)
      gammak = gamma[kk]
      M = diag(1,p)-Pk
      M1 = diag(1,p-1)
      u0 = rep(0,p)
      newM1 = matrix(rep(0,p*(p-1)),p,p-1)
      newM1[-kk,] = M1
      M1 = newM1
      Zpe  = Ak%*%t(X1)
      Wpe = Re(t(M1))%*%t(X1)
      Wk = Re(t(M1))%*%t(X)
    }
    if(model!=1&model!=2)
    {
      ck = C[,kk]
      Pk = ck%*%solve(t(ck)%*%ck)%*%t(ck)
      Ak = solve(t(ck)%*%ck)%*%t(ck)
      Zk = Ak%*%t(X)
      gammak = gamma[kk]
      M = diag(1,p)-Pk
      MM = eigen(M)
      location = which(abs(MM$values)>=10^(-7))
      M1 = MM$vectors[,location]
      Zpe  = Ak%*%t(X1)
      Wpe = Re(t(M1))%*%t(X1)
      Wk = Re(t(M1))%*%t(X)
    }
    time1 = Sys.time()
    ###################################
    modelx = cv.glmnet(t(Wpe),Zpe,nfolds = 5,nlambda=10)
    lambda = modelx$lambda.min
    modelx = glmnet(t(Wpe), Zpe, lambda = lambda)
    pi1 = as.numeric(modelx$beta)
    a1 = modelx$a0
    res = t(Zk) - a1 - t(Wk)%*%pi1
    fisher = mean(res^2)
    res1 = Y - a0 - t(Zk)%*%gammak - t(Wk)%*%Re(t(M1))%*%beta1
    sigmahat = t(res1)%*%res1/n1
    res0 = as.vector(Zk) - modelx$a0 - as.vector(t(Wk)%*%pi1)
    res = Y-a0 - t(Zk)%*%gammak - t(Wk)%*%Re(t(M1))%*%beta1
    U = res*res0
    Un = sum(U)/sqrt(sigmahat*n1*fisher)
    Score_stat = cbind(Score_stat, Un)
    pvalue = 2*(1-pnorm(abs(Un)))
    Pvalue = cbind(Pvalue,pvalue)
    time2 = Sys.time()
#    print(kk)
#    print(c(time1-time0,time2-time1))
  }
  stat1 = max(Score_stat^2)
  size_add = ifelse(stat1>=2*log(r)-log(log(r))-log(3.1415926)-2*log(log(1/0.95)),1,0)#for mcl method
  reject = size_add#for mcl method
  list(stat = Score_stat, reject = reject, maxstat = stat1)
}

################
Infer_decor_2 = function(X,Y,seed,spilit)
{
  X = X
  Y = Y
  n = dim(X)[1]
  p = dim(X)[2]
Omega = diag(1,p)

  if(spilit == "Once")
  {
    n1 = n/2
    X1 = X[1:n1,]
    Y1 = Y[1:n1]   
    X = X[-c(1:n1),]
    Y = Y[-c(1:n1)]  
  }
  if(spilit == "Twice")
  {
    n1 = n/2
    X1 = X[-c(1:n1),]
    Y1 = Y[-c(1:n1)]  
    X = X[c(1:n1),]
    Y = Y[c(1:n1)]  
  }
  Z = solve(t(C)%*%Omega%*%C)%*%t(C)%*%Omega%*%t(X)
  W = t(X)-C%*%Z
  mY = mean(Y)
  mZ = apply(Z,1,mean)
  mW = apply(W,1,mean)
  Z1 = Z
  W1 = W
  for(k in 1:n1)
  {
    if(r!=1)
    {
      Z1[,k] = Z[,k]- mZ
    }
    if(r==1)
    {
      Z1[k] = Z[k]-mZ
    }
    W1[,k] = W1[,k] - mW
  }
  modelx = cv.glmnet(X1,Y1,nfolds=5)
  lambda = modelx$lambda.min
  modelx = glmnet(X1,Y1,lambda = lambda)
  a0 = modelx$a0
  beta1 = as.numeric(modelx$beta)
  res = Y-a0-t(Z1)%*%gamma
  PC = C%*%solve(t(C)%*%C)%*%t(C)
  M = diag(1,p)-PC
  if(model==2)
  {
    M0 = diag(0,r)
    M2 = diag(1,r)
    M1 = rbind(M0,M2)
  }
  if(model!=1&model!=2)
  {
    MM = eigen(M)
    location = which(abs(MM$values)>=10^(-7))
    M1 = eigen(M)$vectors[,location]
  }
  Zpe  = solve(t(C)%*%C)%*%t(C)%*%t(X1)
  Wpe = Re(t(M1))%*%t(X1)
  Zk  = solve(t(C)%*%C)%*%t(C)%*%t(X)
  Wk = Re(t(M1))%*%t(X) 
  Zpe1 = NULL
  Score_stat = Pvalue = NULL
  for(kk in 1:r)
  {
    Zk1 = t(Zpe)[,kk]
    Zk2 = t(Zk)[,kk]
    time0 = Sys.time()
    modelx = cv.glmnet(t(Wpe),Zk1,nfolds=5,nlambda=10)
    time1 = Sys.time()
    time1 - time0
    lambda = modelx$lambda.min
    modelx = glmnet(t(Wpe),Zk1,lambda = lambda)
    pi1 = as.numeric(modelx$beta)
    a1 = modelx$a0
    res0 = as.vector(Zk2) - a1- as.vector(t(Wk)%*%pi1)
    fisher =mean(res0^2)#sum(Zk2^2)/n1-t(pi1)%*%(apply(t(Wk)*as.vector(Zk2),2,mean))
    res1 = Y-a0-t(Zk)%*%gamma - t(Wk)%*%Re(t(M1))%*%beta1
    sigmahat = t(res1)%*%res1/n1
    res = Y-a0 - t(Zk)%*%gamma - t(Wk)%*%Re(t(M1))%*%beta1
    U = res*res0
    Un = sum(U)/sqrt(sigmahat*n1*fisher)
    Score_stat = cbind(Score_stat, Un)
    pvalue = 2*(1-pnorm(abs(Un)))
    Pvalue = cbind(Pvalue,pvalue)
    time1 = Sys.time()
    #print(kk)
  }
  stat1 = max(Score_stat^2)
  size_add = ifelse(stat1>=2*log(r)-log(log(r))-log(3.1415926)-2*log(log(1/0.95)),1,0)#for mcl method
  reject = size_add#for mcl method
  list(stat = Score_stat, reject = reject, maxstat = stat1)
}
###############################
Infer_decor_global = function(X,Y,seed,spilit)
{
  X = X
  Y = Y
  n = dim(X)[1]
  p = dim(X)[2]
Omega = diag(1,p)

  if(spilit == "Once")
  {
    n1 = n/2
    X1 = X[1:n1,]
    Y1 = Y[1:n1]   
    X = X[-c(1:n1),]
    Y = Y[-c(1:n1)]  
  }
  if(spilit == "Twice")
  {
    n1 = n/2
    X1 = X[-c(1:n1),]
    Y1 = Y[-c(1:n1)]  
    X = X[c(1:n1),]
    Y = Y[c(1:n1)]  
  }
  Z = solve(t(C)%*%Omega%*%C)%*%t(C)%*%Omega%*%t(X)
  W = t(X)-C%*%Z
  mY = mean(Y)
  mZ = apply(Z,1,mean)
  mW = apply(W,1,mean)
  Z1 = Z
  W1 = W
  for(k in 1:n1)
  {
    if(r!=1)
    {
      Z1[,k] = Z[,k]- mZ
    }
    if(r==1)
    {
      Z1[k] = Z[k]-mZ
    }
    W1[,k] = W1[,k] - mW
  }
  Zk  = solve(t(C)%*%C)%*%t(C)%*%t(X)
  Zpe1 = NULL
  Score_stat = Pvalue = NULL
  for(kk in 1:r)
  {
    Zk2 = t(Zk)[,kk]
    res0 = as.vector(Zk2)
    res = Y-mY
    U = res*res0
    fisher =sum(Zk2^2)/n1
    sigmahat = t(res)%*%res/n1
    Un = sum(U)/sqrt(sigmahat*n1*fisher)
    Score_stat = cbind(Score_stat, Un)
    pvalue = 2*(1-pnorm(abs(Un)))
    Pvalue = cbind(Pvalue,pvalue)
    time1 = Sys.time()
    #print(kk)
  }
  stat1 = max(Score_stat^2)
  size_add = ifelse(stat1>=2*log(r)-log(log(r))-log(3.1415926)-2*log(log(1/0.95)),1,0)#for mcl method
  reject = size_add#for mcl method
  list(stat = Score_stat, reject = reject, maxstat = stat1)
}
################
Infer_decor_full = function(X,Y)
{
  
  X1 = X
  Y1 = Y

  n = dim(X)[1]
  p = dim(X)[2]
Omega = diag(1,p)

  Z = solve(t(C)%*%Omega%*%C)%*%t(C)%*%Omega%*%t(X)
  W = t(X)-C%*%Z
  mY = mean(Y)
  mZ = apply(Z,1,mean)
  mW = apply(W,1,mean)
  Z1 = Z
  W1 = W
  n1 = n
  for(k in 1:n1)
  {
    if(r!=1)
    {
      Z1[,k] = Z[,k]- mZ
    }
    if(r==1)
    {
      Z1[k] = Z[k]-mZ
    }
    W1[,k] = W1[,k] - mW
  }
  res = Y-mY-t(Z1)%*%gamma
  modelx = cv.glmnet(X,Y,nfolds=5)
  lambda = modelx$lambda.min
  modelx = glmnet(X,Y,lambda = lambda)
  a0 = modelx$a0
  beta1 = as.numeric(modelx$beta)
  res = Y-a0-t(Z1)%*%gamma - t(W1)%*%beta1 
  Score_stat = Pvalue = NULL
  for(kk in 1:r)
  {
    time0 = Sys.time()
    ck = C[,kk]
    Pk = ck%*%solve(t(ck)%*%ck)%*%t(ck)
    Ak = solve(t(ck)%*%ck)%*%t(ck)
    Zk = Ak%*%t(X)
    gammak = gamma[kk]
    M = diag(1,p)-Pk
    ##########################
    if(model==1|model==2)
    {
      M1 = diag(1,p-1)
      u0 = rep(0,p)
      newM1 = matrix(rep(0,p*(p-1)),p,p-1)
      newM1[-kk,] = M1
      M1 = newM1
    }
    if(model!=1&model!=2)
    {
      ck = C[,kk]
      Pk = ck%*%solve(t(ck)%*%ck)%*%t(ck)
      Ak = solve(t(ck)%*%ck)%*%t(ck)
      Zk = Ak%*%t(X)
      gammak = gamma[kk]
      M = diag(1,p)-Pk
      MM = eigen(M)
      location = which(abs(MM$values)>=10^(-7))
      M1 = MM$vectors[,location]
    }
    ###################################
    Zpe  = Ak%*%t(X)
    Wpe = Re(t(M1))%*%t(X)
    Wk = Re(t(M1))%*%t(X)
    modelx = cv.glmnet(t(Wpe),Zpe,nfolds = 5,nlambda=10)
    lambda = modelx$lambda.min
    modelx = glmnet(t(Wpe), Zpe, lambda = lambda)
    pi1 = as.numeric(modelx$beta)
    a1 = modelx$a0
    res0 = as.vector(Zk) - a1 -  as.vector(t(Wk)%*%pi1)
    fisher =mean(res0^2)#sum(Zpe^2)/n1-t(pi1)%*%(apply(t(Wpe)*as.vector(Zpe),2,mean))
    res1 = Y - a0-t(Wpe)%*%Re(t(M1))%*%beta1
    sigmahat = t(res1)%*%res1/n1
    res = Y-a0 - t(Zk)%*%gammak - t(Wk)%*%Re(t(M1))%*%beta1
    U = res*res0
    Un = sum(U)/sqrt(sigmahat*n1*fisher)
    Score_stat = cbind(Score_stat, Un)
    pvalue = 2*(1-pnorm(abs(Un)))
    Pvalue = cbind(Pvalue,pvalue)
    time1 = Sys.time()
    print(kk)
    print(time1-time0)
  }
  stat1 = max(Score_stat^2)
  size_add = ifelse(stat1>=2*log(r)-log(log(r))-log(3.1415926)-2*log(log(1/0.95)),1,0)#for mcl method
  reject = size_add#for mcl method
  list(stat = Score_stat, reject = reject, maxstat = stat1)
}
#################
Infer_decor_full_2 = function(X,Y)
{
  
  X = X
  Y = Y

  n = dim(X)[1]
  p = dim(X)[2]
Omega = diag(1,p)

  n1 = n
  Z = solve(t(C)%*%Omega%*%C)%*%t(C)%*%Omega%*%t(X)
  W = t(X)-C%*%Z
  mY = mean(Y)
  mZ = apply(Z,1,mean)
  mW = apply(W,1,mean)
  Z1 = Z
  W1 = W
  for(k in 1:n1)
  {
    if(r!=1)
    {
      Z1[,k] = Z[,k]- mZ
    }
    if(r==1)
    {
      Z1[k] = Z[k]-mZ
    }
    W1[,k] = W1[,k] - mW
  }
  modelx = cv.glmnet(X,Y,nfolds=5)
  lambda = modelx$lambda.min
  modelx = glmnet(X,Y,lambda = lambda)
  a0 = modelx$a0
  beta1 = as.numeric(modelx$beta)
  res = Y-a0-t(Z)%*%gamma
  PC = C%*%solve(t(C)%*%C)%*%t(C)
  M = diag(1,p)-PC
  if(model==2)
  {
    M0 = diag(0,r)
    M2 = diag(1,r)
    M1 = rbind(M0,M2)
  }
  if(model!=1&model!=2)
  {
    MM = eigen(M)
    location = which(abs(MM$values)>=10^(-8))
    M1 = eigen(M)$vectors[,location]
  }
  Zpe  = solve(t(C)%*%C)%*%t(C)%*%t(X)
  Wpe = Re(t(M1))%*%t(X)
  Zk  = solve(t(C)%*%C)%*%t(C)%*%t(X)
  Wk = Re(t(M1))%*%t(X) 
  Zpe1 = NULL
  Score_stat = Pvalue = NULL
  for(kk in 1:r)
  {
    Zk1 = t(Zpe)[,kk]
    Zk2 = t(Zk)[,kk]
    modelx = cv.glmnet(t(Wpe),Zk1,nfolds=5)
    lambda = modelx$lambda.min
    time0 = Sys.time()
    modelx = glmnet(t(Wpe),Zk1,lambda = lambda,nlambda=10)
    time1 = Sys.time()
    time1-time0
    pi1 = as.numeric(modelx$beta)
    a1 = modelx$a0
    res0 = as.vector(Zk2) - a1 - as.vector(t(Wk)%*%pi1)
    fisher = mean(res0^2)#sum(Zk1^2)/n1-t(pi1)%*%(apply(t(Wpe)*as.vector(Zk1),2,mean))
    res1 = Y - a0 - t(Wpe)%*%Re(t(M1))%*%beta1
    sigmahat = t(res1)%*%res1/n1
    res = Y-a0 - t(Zk)%*%gamma - t(Wk)%*%Re(t(M1))%*%beta1
    U = res*res0
    Un = sum(U)/sqrt(sigmahat*n1*fisher)
    Score_stat = cbind(Score_stat, Un)
    pvalue = 2*(1-pnorm(abs(Un)))
    Pvalue = cbind(Pvalue,pvalue)
    time1 = Sys.time()
    #print(kk)
  }
  stat1 = max(Score_stat^2)
  size_add = ifelse(stat1>=2*log(r)-log(log(r))-log(3.1415926)-2*log(log(1/0.95)),1,0)#for mcl method
  reject = size_add#for mcl method
  list(stat = Score_stat, reject = reject, maxstat = stat1)
}

###################
Infer_decor_full_global = function(X,Y)
{
  
  X = X
  Y = Y
  n = dim(X)[1]
  p = dim(X)[2]
Omega = diag(1,p)

  n1 = n
  Z = solve(t(C)%*%Omega%*%C)%*%t(C)%*%Omega%*%t(X)
  W = t(X)-C%*%Z
  mY = mean(Y)
  mZ = apply(Z,1,mean)
  mW = apply(W,1,mean)
  Z1 = Z
  W1 = W
  for(k in 1:n1)
  {
    if(r!=1)
    {
      Z1[,k] = Z[,k]- mZ
    }
    if(r==1)
    {
      Z1[k] = Z[k]-mZ
    }
    W1[,k] = W1[,k] - mW
  }
  Zk  = solve(t(C)%*%C)%*%t(C)%*%t(X)
  Zpe1 = NULL
  Score_stat = Pvalue = NULL
  for(kk in 1:r)
  {
    Zk2 = t(Zk)[,kk]
    res0 = as.vector(Zk2)
    res = Y-mY-t(Zk)%*%gamma
    fisher = sum(res0^2)/n1
    sigmahat = sum(res^2)/n1
    U = res*res0
    Un = sum(U)/sqrt(sigmahat*n1*fisher)
    Score_stat = cbind(Score_stat, Un)
    pvalue = 2*(1-pnorm(abs(Un)))
    Pvalue = cbind(Pvalue,pvalue)
    time1 = Sys.time()
    #print(kk)
  }
  stat1 = max(Score_stat^2)
  size_add = ifelse(stat1>=2*log(r)-log(log(r))-log(3.1415926)-2*log(log(1/0.95)),1,0)#for mcl method
  reject = size_add#for mcl method
  list(stat = Score_stat, reject = reject, maxstat = stat1)
}