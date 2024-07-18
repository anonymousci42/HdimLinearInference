source("Functions.R")#set the working directory as the path of this script
###############Experiment Setting
##n: sample size
##p: dimension of covariates
##M: Working matrix. "I"-Identity matrix; "Hat"-Estimated precision matrix; "True"-True precision matrix Omega.
##deltas: the sequence of the value in the alternative
##covstru:Sigma. "AR" AR(1) structure with rho. "BCS"-Blockwise compound symmetry with rho.
##method: construction method of J_{PE}. 1: first method; 2: second method.
##combine: "False": single data splitting; "True": double data splitting.
nn=n= 200
p= 600
model= 1
size_pe = size_mcl = 0
deltas= 0
covstru="AR"
M = "I"
method= 1
combine="False"
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
record_stat_tpe = record_stat_npe = Tstat = PE_stat = Stats_all =  NULL
for(i in 1:500)
{
  ##############Continuous Datadelta
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
  if(model == 3)
  {
    r = p/2 + 2
    beta = rep(0,p)
    gamma = rep(0,r)
    C = matrix(rep(0,r*p),p,r)
    for(K in 1:(p/2))
    {
      C[K,K] = 1
    }
    C[p-1,r-1] = 1
    #  beta[p-1] = 1
    beta[p] = delta
    C[p,r-1] = -1
    C[p-1,r] = 1
    #  beta[p-1] = 1
    #  beta[p-2] = 1
    C[p-2,r] = -1
    t(C)%*%beta
  }
  if(model == 4)
  {
    r = p/2
    beta = rep(0,p)
    index = 6*c(1:15)
    beta[index] =delta*rep(1,15)
    gamma = rep(0,r)
    C = matrix(rep(0,r*p),p,r)
    for(K in 1:r)
    {
      C[K,K] = 1
      C[K+p/2,K] = -1
    }
  }
    if(model == 6)
    {
      r = p/2
      beta = rep(0,p)
      index = 6*c(1:15)
      beta[index] =delta*rep(1,15)
      gamma = rep(0,r)
      set.seed(1)
      C = matrix(rnorm(r*p,0,1),p,r)/sqrt(r)
      set.seed(NULL)
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
  #######################Discrete data case
  if(model == 5)
  {
    seeds = 1234+i
    kk = 300
    X1 = matrix(rnorm(kk*n,0,1),n,kk)
    X1 = ifelse(X1<=qnorm(1/3,0,1),1,0) + ifelse(X1>qnorm(2/3,0,1),2,0)
    p1 = p - kk
    X2 = 1*gen_error(n,p1,rho=0.5)
    X = cbind(X1,X2)
    index = 1:2
    beta = rep(0,p)
    beta[index] = 0.5
    beta[(kk+1):(kk+2)] = 1
    epsilon = rnorm(nn,0,1)
    Y = X%*%beta + epsilon
    Y = as.vector(Y)
    ###################
    Z1 = ifelse(X1==1,1,0)
    Z2 = ifelse(X1==2,1,0)
    Xnew = cbind(Z1,Z2,X2)
    p2 = dim(Xnew)[2]
    n = nn = dim(Xnew)[1]
    ###################
    r = kk
    gamma = rep(0,r)
    gamma[index] = delta
    C = matrix(rep(0,r*p2),p2,r)
    for(K in 1:kk)
    {
      C[K,K] = 2
      C[K+r,K] = -1
    }
    X = Xnew
    data_X = Xnew
    data_Y = Y
    Omega = diag(1,p)
  }
  ########################
  if(combine!="True")
  {
    T_PE1 =T_PE2 =T_PE3 = Tn1 =Tn2 = Tn3 = 0
    Result1 = Inference_Npe(data_X,data_Y,CovMatrix,"Once",M,1234+i)
    npe1 = Result1$stat_npe
    Tn1 = npe1
    if(method == 1)
    {
      Result = Infer_decor_1(data_X,data_Y,seeds,"Once")
    }
    if(method == 2&model!=1)
    {
      Result = Infer_decor_2(data_X,data_Y,seeds,"Once")
    }
    if(method == 2&model==1)
    {
      Result = Infer_decor_global(data_X,data_Y,seeds,"Once")
    }
    maxstat = Result$maxstat
    PE = Result$stat
    size_mcl = size_mcl + Result$reject
    cutoff = 2*log(r) + 2*sqrt(log(n/2))*log(log(r))
    T_PE1 = abs(Tn1) + sqrt(r)*sum(ifelse(PE^2>=cutoff,PE^2,0))
    size_pe = size_pe + ifelse(maxstat>=cutoff,1,0)
  }
  if(combine == "True")
  {
    T_PE1 =T_PE2 =T_PE3 = Tn1 =Tn2 = Tn3 = 0
    Result1 = Inference_Npe(data_X,data_Y,CovMatrix,cv="Once",M,1234+i)
    Result2 = Inference_Npe(data_X,data_Y,CovMatrix,cv="Twice",M,1234+i)
    npe1 = Result1$stat_npe
    npe2 = Result2$stat_npe
    Tn2 = 1/sqrt(2)*(npe1 + npe2)
    Result3 = Inference_Npe(data_X,data_Y,CovMatrix,cv="Full",M,1234+i)
    Tn3 = Result3$stat_npe
    if(method == 1)
    {
      Result = Infer_decor_full(data_X,data_Y)
    }
    if(method == 2&model!=1)
    {
      Result = Infer_decor_full_2(data_X,data_Y)
    }
    if(method == 2&model==1)
    {
      Result = Infer_decor_full_global(data_X,data_Y)
    }
    maxstat = Result$maxstat
    PE = Result$stat
    size_mcl = size_mcl + Result$reject
    cutoff = 2*log(r) + 2*sqrt(log(n))*log(log(r))
    T_PE2 = abs(Tn2) + sqrt(r)*sum(ifelse(PE^2>=cutoff,PE^2,0))
    T_PE3 = abs(Tn3) + sqrt(r)*sum(ifelse(PE^2>=cutoff,PE^2,0))
  }
  Tstat = c(Tstat,maxstat)
  Stat_add = c(Tn1,Tn2,Tn3,T_PE1,T_PE2,T_PE3)
  Stats_all = rbind(Stats_all, Stat_add)
  size_all = apply(ifelse(abs(Stats_all)>=qnorm(0.975),1,0),2,mean)
  PE_stat = rbind(PE_stat,PE)
  size1 = size_mcl/i
  information = c(i,maxstat, size_all, size1)
  print(information)
  record_stat = data.frame(Maxstat=Tstat)
  write.csv(information,paste("Size_method_",method,"_covmatrix_",covstru,"_combine_",combine,"_n_",n,"_p_",p,"_model_",model,"_delta_",delta*10,"_.csv",sep=""),row.names = FALSE,col.names = FALSE)
 # write.csv(PE_stat,paste("PEstat_method_",method,"_covmatrix_",covstru,"_combine_",combine,"_n_",n,"_p_",p,"_model_",model,"_delta_",delta*10,"_.csv",sep=""),row.names = FALSE,col.names = FALSE)
 # write.csv(record_stat,paste("Maxstat_method_",method,"_covmatrix_",covstru,"_combine_",combine,"_n_",n,"_p_",p,"_model_",model,"_delta_",delta*10,"_.csv",sep=""),row.names = FALSE,col.names = FALSE)
 # write.csv(Stats_all,paste("Basicstat_method_",method,"_covmatrix_",covstru,"_combine_",combine,"_n_",n,"_p_",p,"_model_",model,"_delta_",delta*10,"_.csv",sep=""),row.names = FALSE,col.names = FALSE)
}
}
#}