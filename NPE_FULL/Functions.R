library(glmnet)
library(ncvreg)
library(MASS)
library(SIHR)
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
scad1 = function(lambda,beta)
{
  a = 3.7
  scaddot = ifelse(beta<=lambda,lambda,(max(a*lambda-abs(beta),0))/(a-1))
  return(scaddot)
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

#########################Partial Penalized Method of Shi. et.al (2019)

Infer_partial_pen = function(data_X,data_Y,C0)
{
  active = 1:(dim(C0)[1])
  #C = diag(length(active))
  C = t(C0)
  X = data_X
  Y = data_Y
  n = dim(X)[1]
  p = dim(X)[2]
  M = rep(1,p)==0
  M[active] = 1
  pen_factor = rep(1,p)
  pen_factor[M!=0] = 0
  #  tmp = cv.ncvreg(data_X,data_Y,nfolds=5,penalty = "SCAD",nlambda = 50,penalty.factor = pen_factor)
  modelx = cv.glmnet(data_X,data_Y,nfolds=5,nlambda = 50,penalty.factor = pen_factor)
  lambda = modelx$lambda.min
  modelx = glmnet(data_X,data_Y,lambda = lambda,penalty.factor = pen_factor)
  beta_ini = as.numeric(modelx$beta)
  ###########One step sparse estimator of scad
  pen_factor_1 = scad1(lambda,beta_ini)
  pen_factor_1[M!=0] = 0
  #  modelx = cv.glmnet(data_X,data_Y,nfolds=5,nlambda = 50,penalty.factor = pen_factor)
  #  lambda = modelx$lambda.min
  modelx = glmnet(data_X,data_Y,lambda = lambda,penalty.factor = pen_factor_1)
  beta_ini = as.numeric(modelx$beta)
  a0 = modelx$a0
  # beta_ini <- tmp$fit$beta[,tmp$min]
  S = which(beta_ini!=0)
  MM = which(M!=0)
  rindex = c(MM,S[!S%in%(MM)])
  mm = length(MM)
  rm = length(rindex)
  resid = as.vector(Y - X%*%beta_ini - a0)
  sig = mean(resid^2)
  Z = X[,rindex]
  Zr = Z
  D0 = t(Z)%*%(Z)/n
  D3 = sig*t(Zr)%*%Zr/n
  Acov = solve(D0)%*%D3%*%solve(t(D0))
  Acovm = Acov[1:mm,1:mm]
  ########Statistic construction
  betahat = beta_ini#[-1]
  betahat = betahat[active]
  a1 = C%*%betahat - gamma
  a2 = (C)%*%Acovm%*%t(C)
  Tn_ini = n*t(a1)%*%solve(a2)%*%a1#sqrt(n*T)*betahat[2]/sqrt(Acovm[2,2])#
  rej = ifelse(abs(Tn_ini)>=qchisq(0.95,dim(C0)[2]),1,0)
  list(Tn_ini = Tn_ini,rej=rej)
}