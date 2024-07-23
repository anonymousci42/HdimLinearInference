source("Functions.R")
###############Experiment Setting
##n: sample size
##p: dimension of covariates
##M: Working matrix. "I"-Identity matrix; "Hat"-Estimated precision matrix; "True"-True precision matrix Omega.
##deltas: the sequence of the value in the alternative
##covstru:Sigma. "AR" AR(1) structure with rho. "BCS"-Blockwise compound symmetry with rho.
##combine: "False" single data splitting; "True": double data splitting
nn=n= 200
p= 600
M = "I"
model= 5#For the comparison of PPE, only use model n5
size_pe = size_mcl = 0
deltas=seq(0,1,length=6)
covstru="AR"
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
  record_stat_tpe = record_stat_npe = Tstat = PE_stat = Stats_all = Tstat_pp = NULL
for(i in 1:500)
{
  if(model==5)
  {
    seeds = 1234+i
    kk = 0.4*n
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
    Omega = diag(1,p)
    
  }
  #####################
  ########################
  if(combine!="True")
  {
    Tn1 =Tn2 = Tn3 = 0
    Result1 = Inference_Npe(data_X,data_Y,CovMatrix,"Once",M,1234+i)
    npe1 = Result1$stat_npe
    Tn1 = npe1
  }
  if(combine == "True")
  {
    Tn1 =Tn2 = Tn3 = 0
    Result1 = Inference_Npe(data_X,data_Y,CovMatrix,cv="Once",M,1234+i)
    Result2 = Inference_Npe(data_X,data_Y,CovMatrix,cv="Twice",M,1234+i)
    npe1 = Result1$stat_npe
    npe2 = Result2$stat_npe
    Tn2 = 1/sqrt(2)*(npe1 + npe2)
    Result3 = Inference_Npe(data_X,data_Y,CovMatrix,cv="Full",M,1234+i)
    Tn3 = Result3$stat_npe
  }
  Stat_add = c(Tn1,Tn2,Tn3)
  Stats_all = rbind(Stats_all, Stat_add)
  size_all = apply(ifelse(abs(Stats_all)>=qnorm(0.975),1,0),2,mean)
  #####################Comparison
  Result_pp = Infer_partial_pen(data_X,data_Y,C = C0)
  Tstat_pp = c(Tstat_pp, Result_pp$Tn_ini)
  size_pp = mean(ifelse(Tstat_pp>=qchisq(0.95,df = r),1,0))
  #########################################
  information = c(i,size_all,size_pp)
  print(information)
  write.csv(information,paste("Size","_covmatrix_",covstru,"_combine_",combine,"_n_",n,"_p_",p,"_model_",model,"_delta_",delta*10,"_.csv",sep=""),row.names = FALSE,col.names = FALSE)
 # write.csv(PE_stat,paste("PEstat_method_",method,"_covmatrix_",covstru,"_combine_",combine,"_n_",n,"_p_",p,"_model_",model,"_delta_",delta*10,"_.csv",sep=""),row.names = FALSE,col.names = FALSE)
 # write.csv(record_stat,paste("Maxstat_method_",method,"_covmatrix_",covstru,"_combine_",combine,"_n_",n,"_p_",p,"_model_",model,"_delta_",delta*10,"_.csv",sep=""),row.names = FALSE,col.names = FALSE)
 # write.csv(Stats_all,paste("Basicstat","_covmatrix_",covstru,"_combine_",combine,"_n_",n,"_p_",p,"_model_",model,"_delta_",delta*10,"_.csv",sep=""),row.names = FALSE,col.names = FALSE)
 # write.csv(Tstat_pp,paste("Stat_PP","_covmatrix_",covstru,"_combine_",combine,"_n_",n,"_p_",p,"_model_",model,"_delta_",delta*10,"_.csv",sep=""),row.names = FALSE,col.names = FALSE)
 }
}
