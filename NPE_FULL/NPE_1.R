source("Functions.R")#set the working directory as the path of this script
###############Experiment Setting
##n: sample size
##p: dimension of covariates
##M: Working matrix. "I"-Identity matrix; "Hat"-Estimated precision matrix; "True"-True precision matrix Omega.
##deltas: the sequence of the value in the alternative
##covstru:Sigma. "AR" AR(1) structure with rho. "BCS"-Blockwise compound symmetry with rho.
nn=n= 200
p= 600
M="I"
model= 1
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
record_stat_tpe = record_stat_npe = Tstat = Times = NULL
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
    C[p-2,r] = -1
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
  ######################data generating
  nn = n
  covsigma = 1
  X = mvrnorm(n,rep(0,p),covsigma*CovMatrix)
  epsilon = rnorm(nn,0,1)
  Y = X%*%beta + epsilon
  Y = as.vector(Y)
  data_X = X
  data_Y = Y
  Omega = diag(1,p)
  }
  #######################Only for Model N5
  if(model == 5)
  {
    kk = 300
    X1 = matrix(rnorm(kk*n,0,1),n,kk)
    X1 = ifelse(X1<=qnorm(1/3,0,1),1,0) + ifelse(X1>qnorm(2/3,0,1),2,0)
    p1 = p - kk
    X2 = 1*gen_error(n,p1,rho=0.5)
    X = cbind(X1,X2)
    index = 1:2
    beta1 = rep(0,kk)
    beta2 = rep(0,kk)
    beta3 = rep(0,p1)
    beta1[index] = 0.5
    beta2[index] = 2*beta1[index] - delta
    beta3[index] = 1
    beta = c(beta1,beta2,beta3)
    #beta[index] = 0.5
    #beta[(kk+1):(kk+2)] = 1
    epsilon = rnorm(nn,0,1)
    ###################
    Z1 = ifelse(X1==1,1,0)
    Z2 = ifelse(X1==2,1,0)
    Xnew = cbind(Z1,Z2,X2)
    Y = Xnew%*%beta + epsilon
    Y = as.vector(Y)
    p2 = dim(Xnew)[2]
    n = nn = dim(Xnew)[1]
    ###################
    r = kk
    gamma = rep(0,r)
    #gamma[index] = delta
    C = matrix(rep(0,r*p2),p2,r)
    for(K in 1:kk)
    {
      C[K,K] = 2
      C[K+r,K] = -1
    }
    C0 = C[1:(2*kk),]
    X = Xnew
    data_X = Xnew
    data_Y = Y
  }
  ########################
  op <- options(digits.secs = 6)
  ######################single data splitting
  time0 = Sys.time()
  Result1 = Inference_Npe(data_X,data_Y,CovMatrix,"Once",M,1234+i)
  npe1 = Result1$stat_npe
  Tn1 = npe1 #Statistic of single data splitting
  time1 = Sys.time()
  #######################Double data splitting
  Result2 = Inference_Npe(data_X,data_Y,CovMatrix,cv="Twice",M,1234+i)
  npe2 = Result2$stat_npe
  Tn2 = 1/sqrt(2)*(npe1 + npe2)#Statistic of double data splitting
  time2 = Sys.time()
  #######################Full data
  Result3 = Inference_Npe(data_X,data_Y,CovMatrix,cv="Full",M,1234+i)
  Tn3 = Result3$stat_npe#Statistic of full data
  time3 = Sys.time()
  #######################Result summary
  Times = rbind(Times, as.numeric(c(difftime(time1,time0,units="secs"), difftime(time2,time0,units="secs"), difftime(time3,time2,units="secs"))))#Computation time record
  Ts = c(Tn1,Tn2,Tn3)
  Tstat = rbind(Tstat,Ts)
  size_npe = apply(ifelse(abs(Tstat)>=qnorm(0.975),1,0),2,mean)#Empirical rejection rate
  information = c(i,Tn1,size_npe,apply(Times,2,mean))
  print(information)
  record_stat = data.frame(npe=Tstat)
  ####Save the empirical rejection rate of the finished simulations.
  write.csv(information,paste("Size_M_",M,"_covmatrix_",covstru,"_n_",n,"_p_",p,"_model_",model,"_delta_",delta*10,"_.csv",sep=""),row.names = FALSE,col.names = FALSE) 
  ####Save the statistic value of the finished simulations
  write.csv(record_stat,paste("stat_M_",M,"_covmatrix_",covstru,"_n_",n,"_p_",p,"_model_",model,"_delta_",delta*10,"_.csv",sep=""),row.names = FALSE,col.names = FALSE)
}
}
