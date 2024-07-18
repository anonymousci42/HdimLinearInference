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
nullpower_append <- function(dataset, power){
  iter = nrow(dataset)
  dataset = dataset[,-c(1,2)]
  take = ncol(power)
  dataset = dataset[,c(1:(take-1))]
  nnrow = colSums(abs(dataset) > 1.96 * 1)/iter
  nptest = rbind(c(0,nnrow),power)
  return(nptest)
}

z_intercept <- function(x, y, b0 = 0){
  n = ncol(x)
  p = nrow(x)
  beta0 = rep(b0, p)
  z_mat = matrix(data = NA, nrow = p, ncol = (n*(n-1)/2))
  index = Rfast::comb_n(n, 2)
  
  for(i in 1:ncol(index)){
    z_mat[,i] = (x[,index[1,i]] - x[,index[2,i]]) * (y[index[1,i]] - y[index[2,i]] - sum(beta0 *  (x[,index[1,i]] - x[,index[2,i]])))
  }
  
  return(z_mat)
}
#################User functions
Inference2 = function(X,Y,spilit,seednum)#For n=200,p<=400 moderate dimensional cases. Otherwise can not allocate enough memory.
{
  X = X
  Y = as.vector(Y)
  n = dim(X)[1]
  p = dim(X)[2]
  Omega = diag(1,p)#solve(CovMatrix)
  indexX = 1:p
  X1 = X
  Y1 = Y
  n1 = n
  if(spilit == "Once")
  {
    n1 = n/2
    X1 = X[1:n1,]
    Y1 = Y[1:n1]   
  }
  if(spilit == "Twice")
  {
    n1 = n/2
    X1 = X[-c(1:n1),]
    Y1 = Y[-c(1:n1)]   
    
  }
  if(spilit == "Once")
  {
    X = X[-c(1:n1),]
    Y = Y[-c(1:n1)]  
  }
  if(spilit == "Twice")
  {
    X = X[c(1:n1),]
    Y = Y[c(1:n1)]  
  }
  Z = solve(t(C)%*%C)%*%t(C)%*%t(X)
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
  modelx = cv.glmnet(X1,Y1,noflds=5)
  lambda = modelx$lambda.min
  modelx = glmnet(X1,Y1,lambda = lambda)
  beta1 = as.numeric(modelx$beta)
  res = Y-mY-t(Z1)%*%gamma - t(W1)%*%beta1 
  res_oracle = Y-mY-t(Z1)%*%t(C)%*%beta1 - t(W1)%*%beta1
  D = sapply(1:ncol(Z),function(x) Z1[,x] * res[x])
  D_oracle = sapply(1:ncol(Z),function(x) Z1[,x] * res[x])
  if(r==1)
  {
    D = t(matrix(D))  
  }
  size = -1
  apply(X,1,sd)
  ###################
  D_new = D
  n1=n=dim(D)[2]
  #####Power enhancement
  PE1 = 0
  ###########
  record_matrix = rep(list(1),n*(n-1))
  Dstar = rep(0,dim(D)[1])
  index_D = 1:n
  index_r = 1:r
  pointer = pointer1 =  1
  dim(D)
  Tn = 0
  r_index = 1
  trace_est = 0
  for(j in 1:n)
  {
    for(k in index_D[-j])
    { 
      temp1 = t(D[,j])%*%(D[,k])
      temp2 = t(D[,j])%*%(D[,k])
     # gc()
      # record_matrix[[pointer]] = temp1
      Tn  = Tn + as.numeric(temp1)
      # pointer = pointer + 1
      #  D_deduce = D[,-c(j,k)]
      #  Dmean1 = apply(D_deduce,1,mean)
      #  A1 = D[,j]-Dmean1
      #  A2 = D[,k]-Dmean1
      #  trace_temp = A1%*%t(D[,j])%*%A2%*%t(D[,k])
      trace_est = trace_est + as.numeric(temp2^2)
    }
  }
  W = Tn/n#(n*(n-1))
  estTB2 = 2*(trace_est)/(n*(n-1))
  T = W/sqrt(estTB2)
  T_pe = T + PE1 
  size_npe =  ifelse(abs(T)>=qnorm(0.975),1,0)
  size_pe =  ifelse(abs(T_pe)>=qnorm(0.975),1,0)
  list(size_npe = size_npe,size_pe = size_pe,stat_npe = T,stat_pe = T_pe, sigma = trace_est, PE  =PE1)
}
############

Infer_decor_1 = function(X,Y,seed,spilit)
{
  X = X
  Y = Y
  n = dim(X)[1]
  p = dim(X)[2]
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
  modelx = cv.glmnet(X1,Y1,noflds=5)
  lambda = modelx$lambda.min
  modelx = glmnet(X1,Y1,lambda = lambda)
  beta1 = as.numeric(modelx$beta)
  res = Y-mY-t(Z1)%*%gamma - t(W1)%*%beta1 
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
    M1 = eigen(M)$vectors[,c((p-1):(1))]
    #t(M1)%*%M1
    #M1%*%t(M1)
    Zpe  = Ak%*%t(X1)
    Wpe = Re(t(M1))%*%t(X1)
    #t(X1)[,1]
    #Wpe[,1]
    Wk = Re(t(M1))%*%t(X)
    modelx = cv.glmnet(t(Wpe),Zpe,nfolds = 5)
    lambda = modelx$lambda.min
    modelx = glmnet(t(Wpe), Zpe, lambda = lambda)
    pi1 = as.numeric(modelx$beta)
    fisher =sum(Zk^2)/n1-t(pi1)%*%(apply(t(Wk)*as.vector(Zk),2,mean))
    res1 = Y - t(Zk)%*%gammak - t(Wk)%*%Re(t(M1))%*%beta1
    sigmahat = t(res1)%*%res1/n1
    res0 = as.vector(Zk) - as.vector(t(Wk)%*%pi1)
    res = Y-t(Zk)%*%gammak - t(Wk)%*%Re(t(M1))%*%beta1
    U = res*res0
    Un = sum(U)/sqrt(sigmahat*n1*fisher)
    #(X)[1,]
    #recovery = t(Wk)%*%Re(t(M1))
    #recovery[1,]
    Score_stat = cbind(Score_stat, Un)
    pvalue = 2*(1-pnorm(abs(Un)))
    Pvalue = cbind(Pvalue,pvalue)
    time1 = Sys.time()
    #print(kk)
    #  print(time1-time0)
  }
  asd1 = ifelse(Score_stat^2>=2*log(r) + 2*log(log(r)),Score_stat^2,0)
  Strong = sum(asd1)
  num = length(asd1[asd1!=0])
  PE1 =  sqrt(r)*sum(asd1)
  stat1 = max(Score_stat^2)
  size_add = ifelse(stat1>=2*log(r)-log(log(r))-log(3.1415926)-2*log(log(1/0.95)),1,0)
  reject = size_add
  list(Strong = Strong, num = num, PE1 = PE1, reject = reject, maxstat = stat1)
}

################
Infer_decor_2 = function(X,Y,seed,spilit)
{
  X = X
  Y = Y
  n = dim(X)[1]
  p = dim(X)[2]
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
  modelx = cv.glmnet(X1,Y1,noflds=5)
  lambda = modelx$lambda.min
  modelx = glmnet(X1,Y1,lambda = lambda)
  beta1 = as.numeric(modelx$beta)
  res = Y-mY-t(Z1)%*%gamma
  PC = C%*%solve(t(C)%*%C)%*%t(C)
  M = diag(1,p)-PC
  M1 = eigen(M)$vectors[,c((p-r):1)]
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
    modelx = cv.glmnet(t(Wpe),Zk1,nfolds=5)
    lambda = modelx$lambda.min
    modelx = glmnet(t(Wpe),Zk1,lambda = lambda)
    pi1 = as.numeric(modelx$beta)
    fisher =sum(Zk2^2)/n1-t(pi1)%*%(apply(t(Wk)*as.vector(Zk2),2,mean))
    res1 = Y-t(Zk)%*%gamma - t(Wk)%*%Re(t(M1))%*%beta1
    sigmahat = t(res1)%*%res1/n1
    res0 = as.vector(Zk2) - as.vector(t(Wk)%*%pi1)
    res = Y-t(Zk)%*%gamma - t(Wk)%*%Re(t(M1))%*%beta1
    U = res*res0
    Un = sum(U)/sqrt(sigmahat*n1*fisher)
    Score_stat = cbind(Score_stat, Un)
    pvalue = 2*(1-pnorm(abs(Un)))
    Pvalue = cbind(Pvalue,pvalue)
    time1 = Sys.time()
    #print(kk)
  }
  asd1 = ifelse(Score_stat^2>=2*log(r) + 2*log(log(r)),Score_stat^2,0)
  Strong = sum(asd1)
  num = length(asd1[asd1!=0])
  PE1 =  sqrt(r)*sum(asd1)
  stat1 = max(Score_stat^2)
  size_add = ifelse(stat1>=2*log(r)-log(log(r))-log(3.1415926)-2*log(log(1/0.95)),1,0)
  reject = size_add
  list(Strong = Strong, num = num, PE1 = PE1, reject = reject, maxstat = stat1)
}
###############################
Infer_decor_global = function(X,Y,seed,spilit)
{
  X = X
  Y = Y
  n = dim(X)[1]
  p = dim(X)[2]
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
  Wk = Re(t(M1))%*%t(X) 
  Zpe1 = NULL
  Score_stat = Pvalue = NULL
  for(kk in 1:r)
  {
    Zk2 = t(Zk)[,kk]
    res0 = as.vector(Zk2)
    res = Y-t(Zk)%*%gamma
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
  asd1 = ifelse(Score_stat^2>=2*log(r) + 2*log(log(r)),Score_stat^2,0)
  Strong = sum(asd1)
  num = length(asd1[asd1!=0])
  PE1 =  sqrt(r)*sum(asd1)
  stat1 = max(Score_stat^2)
  size_add = ifelse(stat1>=2*log(r)-log(log(r))-log(3.1415926)-2*log(log(1/0.95)),1,0)
  reject = size_add
  list(Strong = Strong, num = num, PE1 = PE1, reject = reject, maxstat = stat1)
}
#######################################
Infer_decor_full_global = function(X,Y)
{
  X = X
  Y = Y
  n = dim(X)[1]
  p = dim(X)[2]
  n1 = n
  Zk  = solve(t(C)%*%C)%*%t(C)%*%t(X)
  Zpe1 = NULL
  Score_stat = Pvalue = NULL
  for(kk in 1:p)
  {
    Zk2 = t(Zk)[,kk]
    res0 = as.vector(Zk2)
    res = Y-t(Zk)%*%gamma
    fisher = sum(res0^2)/n1
    sigmahat = sum(res^2)/n1
    U = res*res0
    Un = sum(U)/sqrt(sigmahat*n1*fisher)
    Score_stat = cbind(Score_stat, Un)
    pvalue = 2*(1-pnorm(abs(Un)))
    Pvalue = cbind(Pvalue,pvalue)
    time1 = Sys.time()
  #  print(kk)
  }
  asd1 = ifelse(Score_stat^2>=2*log(p) + 4*log(log(p)),Score_stat^2,0)
  Strong = sum(asd1)
  num = length(asd1[asd1!=0])
  PE1 =  sqrt(p)*sum(asd1)
  stat1 = max(Score_stat^2)
  size_add = ifelse(stat1>=2*log(p)-log(log(p))-log(3.1415926)-2*log(log(1/0.95)),1,0)
  reject = size_add
  list(Strong = Strong, num = num, PE1 = PE1, reject = reject, maxstat = stat1, stat = Score_stat)
}
# Zhong and Chen 2011 comparison
zhongchen2011 <- function(X, y, beta = 0, delta = 0, Sigma = NULL, small.sig = NULL, T0 = FALSE){
  # elements needed for a statistic
  n = ncol(X)
  p = nrow(X)
  small.sigma = var(y)
  beta0 = c(rep(delta,beta), rep(0, p - beta)) 
  
  ### Calculating the Trace of Sigma.hat^2
  # Gets the matrices that represent x_i^T x_j where i != j
  mlist =  apply(X, FUN = function(x) colSums(x * X), MARGIN = 2)
  diag(mlist) = 0
  
  ## Trace of Sigma.hat^2 and T_np at the same time
  # ijkl cases and storage stuff for trace
  ijk = rep(0, n) # Empty vector used for storage
  ijkl.store = rep(0,n)
  sm = sum(mlist)
  cs = Rfast::colsums(mlist)
  
  # storage for T_np
  diffy = matrix(data = 0, nrow = n, ncol = n)
  dxy = matrix(data = 0, nrow = p, ncol = n)
  dxy2 = matrix(data = 0, nrow = p, ncol = n)
  horizontal_all_arrays_sum = rep(0, p)
  vertical_all_arrays_sum = matrix(data = 0, nrow = p, ncol = n)
  horizontal_within_arrays = matrix(data = 0, nrow = p, ncol = n)
  
  for(i in 1:n){
    ## Trace of Sigma stuff
    ijk[i] = sum(mlist[i,] %*% mlist[,-i])
    ijkl.store[i] = sum(mlist[i,] * sm) - sum(2 * mlist[i,] * cs) - sum(2*mlist[i,] * cs[i]) + sum(2*mlist[i,]^2)
    
    ## T_np stuff
    ## Getting every (X_ia - X_ib) and (y_ia - y_ib - (X_ia - X_ib)^T * Beta0) combination
    dxy2 = X[,i] - X
    diffy[i,] = y[i] - y - Rfast::colsums(dxy2 * beta0)
    
    # one weird fix
    dxy2[ ,1:i] = 0 
    diffy[i,1:i] = 0 
    
    ## Combines the two together
    dxy2 = t(diffy[i,] * t(dxy2))
    
    ## all the storage stuff
    horizontal_within_arrays[,i] = Rfast::rowsums(dxy2)
    vertical_all_arrays_sum = vertical_all_arrays_sum + dxy2
    
    # this replaces the sum(dxy^2) below from when it was an array
    dxy = dxy + dxy2^2
  }
  
  
  # Summations to get the final values for the variance equation
  sijk = sum(ijk)
  sijkl = sum(ijkl.store)
  ij = sum(mlist^2)
  
  tr.Sigma.hat2 = (1/(n * (n-1))) * ij - 
    (2/(n*(n-1)*(n-2))) * sijk + 
    (1/(n*(n-1)*(n-2)*(n-3))) * sijkl
  
  # gets the product of all the (X_ia - X_ib)^T (X_ic - X_id)(y_ia - y_ib - (X_ia - X_ib)^T * Beta0)(y_ic - y_id - (X_ic - X_id)^T * Beta0)
  horizontal_all_arrays_sum = Rfast::rowsums(horizontal_within_arrays)
  
  newresult = sum(horizontal_all_arrays_sum^2) - sum(vertical_all_arrays_sum^2) -  sum(horizontal_within_arrays^2) - 2*sum(vertical_all_arrays_sum *horizontal_within_arrays) + sum(dxy)
  
  # Zhong and Chen statistic
  T = (1/4) * (1/3) * (4*3*2)/(n*(n-1)*(n-2)*(n-3)) * newresult/2 #divide by two because the newresult construction has (x_1 - x_2)(x_3- x_4) and (x_3 - x_4)(x_1 - x_2) both, so you need to halve it
  # varT = (2/(n*(n-1)))*tr.Sigma.hat2 * small.sigma^2
  normalized.T = (n*T)/(sqrt(2*tr.Sigma.hat2) * small.sigma)
  
  if(T0 == TRUE){
    # Calculating the Trace of the Sigma^2 here
    holder = rep(0, nrow(Sigma))
    for(i in 1:nrow(Sigma)){
      holder[i] = sum(Sigma[i,] * Sigma[,i])
    }
    tr.Sigma2 = sum(holder)
    
    # Zhong and Chen statistic with true sigma 
    normalized.T0 = (n*T)/(sqrt(2*tr.Sigma2) * small.sig)
    return(c(T, normalized.T, normalized.T0))
  }
  
  return(c(T,normalized.T))
}
