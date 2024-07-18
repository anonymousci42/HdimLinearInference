source("Functions.R")
###############Experiment Setting
##n: sample size
##p: dimension of covariates
##M: Working matrix of SIHR method, not the same meaning of other script.
##deltas: the sequence of the value in the alternative
##covstru:Sigma. "AR" AR(1) structure with rho. "BCS"-Blockwise compound symmetry with rho.
nn=n= 200
p= 600
M="I"
model= 1 #For the SIHR method, only models 1 and 2 are used.
size_pe = size_npe = 0
deltas=seq(0,1,length=6)
covstru="AR"
###################
if(covstru == "AR")
{
  rho=0.5
  CovMatrix = matrix(rep(1,p^2),p,p)
  index1 = 1:p
  for(i in 1:p)
  {
    for(j in index1[-i])
    {
      CovMatrix[i,j] = rho^(abs(i-j))  
    }
  }
  Omega = solve(CovMatrix)
}

if(covstru == "BCS")
{
  rho=0.7
  covsigma = 1
  p1=p/20
  CovMatrix = matrix(rep(1,p1^2),p1,p1)
  index1 = 1:p1
  for(i in 1:p1)
  {
    for(j in index1[-i])
    {
      CovMatrix[i,j] = rho
    }
  }
  CovMatrix1 = CovMatrix
  for(k in 1:19)
  {
    CovMatrix1 = bdiag(CovMatrix1,CovMatrix)
  }
  CovMatrix1 = as.matrix(CovMatrix1)
  Sigma = CovMatrix1
  CovMatrix = CovMatrix1   
}
############################Numerical Results
for(delta in deltas)
{
  record_stat_tpe = record_stat_npe = Tstat = PE_stat = Stats_all= Rej =NULL
    for(i in 1:500)
{
 ##############Data generation of model N1-N4
  if(model!=5)
  {
  if(model == 1)
  {
    r = p
    beta = rep(0,p)
    index = c(1,2)
    beta[index] = delta*rep(1,2)
    gamma = rep(0,r)
    C = diag(1,p)
  }
  if(model == 2)
  {
    r = p/2
    beta = rep(0,p)
    index = c(1,2,5,r+1,r+5)
    beta[index] = delta*rep(1,5)
    gamma = rep(0,r)
    C = matrix(rep(0,r*p),p,r)
    for(K in 1:r)
    {
      C[K,K] = 1
    }
  }
  ######################data generating
  seeds = 1234+i
  nn = n
  set.seed(seeds)
  covsigma = 1
  X = mvrnorm(n,rep(0,p),covsigma*CovMatrix)
  epsilon = rnorm(nn,0,1)
  Y = X%*%beta + epsilon
  Y = as.vector(Y)
  data_X = X
  data_Y = Y
  Omega = diag(1,p)
  }
  ########################
  modelx = cv.glmnet(data_X,data_Y)
  beta_ini = as.numeric(coef(modelx,"lambda.min"))
  if(M=="I")
  {
  SIHR = QF(data_X,data_Y,G = 1:r,A = diag(r),model = "linear",beta.init = beta_ini)
  }
  if(M!="I")
  {
    SIHR = QF(data_X,data_Y,G = 1:r,model = "linear",beta.init = beta_ini)
  }
  ## Compute confidence intervals
  CI = ci(SIHR, alpha=0.05, alternative="two.sided")
  Rej = rbind(Rej, ifelse(0!=CI[,2],1,0))
  information = c(i,apply(Rej,2,mean))
  print(information)
  write.csv(apply(Rej,2,mean),paste("Size_SIHR_M_",M,"_covmatrix_",covstru,"_n_",n,"_p_",p,"_model_",model,"_delta_",delta*10,"_.csv",sep=""),row.names = FALSE,col.names = FALSE)
  write.csv(Rej,paste("Rej_SIHR_M_",M,"_covmatrix_",covstru,"_n_",n,"_p_",p,"_model_",model,"_delta_",delta*10,"_.csv",sep=""),row.names = FALSE,col.names = FALSE)
 }
}