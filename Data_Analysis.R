library(tidyr)
library(stringr)
library(ggplot2)
############################################
#set the path of this script firstly.
root_path = getwd()
setwd(paste(root_path,"/PE_FULL/Data_Results_PE",sep=""))
file_names <- dir()
d = length(file_names)
Record  = NULL
for(i in 1:d)
{
  selects = file_names[i]#BS 200 400
  split =  as.vector(str_split(selects,"_"))[[1]]
  method = as.numeric(split[[3]])
  covmatrix = split[5]
  combine = split[7]
  n = as.numeric(split[9])
  p = as.numeric(split[11])
  model = as.numeric(split[13])
  delta = as.numeric(split[15])/10
  data = do.call(rbind,lapply(selects,read.csv)) 
  ERR = round(t(data)[3:9], 3)
  if(method==2|combine=="False"|model==3|model==4|model==5)
  {
    ERR[7]= "NA"
  }
  information = c(model,n,p,delta,combine,covmatrix,method)
  ERR = c(information,ERR)
  Record = rbind(Record,ERR)
}
Record = data.frame(Record)
colnames(Record) = c("model","n","p","delta","combine","covmatrix","method","Tn1","Tn2","Tn3","T_PE1","T_PE2","T_PE3","MCL")
#######################
setwd(paste(root_path,"/PE_FULL/Summary",sep=""))
for(model in 1:5)
{
  for(combine in c("False","True"))
  {
    Record_model = NULL
  for(n in c(200))
    for(covmatrix in c("AR","BCS"))
    for(method in c(1,2))
    for(p in c(600,1200))
    {   
          {
            if(combine =="False")
            {
            Record_sub = Record[Record$n==n&Record$model==model&Record$p==p&Record$combine==combine&Record$covmatrix==covmatrix&Record$method==method,]
            Record_sub1 = Record_sub[order(as.numeric(Record_sub$delta)),-c(8:10,12:14)]
            Record_model = rbind(Record_model,Record_sub1)
            }
          if(combine !="False")
          {
            Record_sub = Record[Record$n==n&Record$model==model&Record$p==p&Record$combine==combine&Record$covmatrix==covmatrix&Record$method==method,]
            Record_sub1 = Record_sub[order(as.numeric(Record_sub$delta)),-c(8:11)]
            Record_model = rbind(Record_model,Record_sub1)
          }
          }
    }
  Record_model = as.data.frame(t(Record_model))
  write.csv(Record_model,paste("ERR_model_",model,"_combine_",combine,".csv",sep=""))
  }
}
#############################For all SIHR Results
setwd(paste(root_path,"/NPE_FULL/SIHR_Results",sep=""))
file_names <- dir()
d = length(file_names)
Record  = NULL
for(i in 1:d)
{
  selects = file_names[i]#BS 200 400
  split =  as.vector(str_split(selects,"_"))[[1]]
  M = split[[4]]
  covmatrix = split[6]
  p = as.numeric(split[10])
  model = as.numeric(split[12])
  delta = as.numeric(split[14])/10
  data = do.call(rbind,lapply(selects,read.csv)) 
  ERR = round(t(data)[3], 3)
  information = c(model,p,delta, M,covmatrix)
  ERR = c(information,ERR)
  Record = rbind(Record,ERR)
}
Record = data.frame(Record)
colnames(Record) = c("model","p","delta","M","covmatrix","SIHR")
########################
setwd(paste(root_path,"/NPE_FULL/Summary",sep=""))
for(model in 1:2)
{    
  for(M in c("I"))
  {
    Record_model = NULL
    for(covmatrix in c("AR","BCS"))
          for(p in c(600,1200))
          {   
            {
                Record_sub = Record[Record$model==model&Record$p==p&Record$M==M&Record$covmatrix==covmatrix,]
                Record_sub1 = Record_sub[order(as.numeric(Record_sub$delta)),]
                Record_model = rbind(Record_model,Record_sub1)
            }
          }
    Record_model = as.data.frame(t(Record_model))
    write.csv(Record_model,paste("SIHR_model_",model,"_M_",M,".csv",sep=""))
  }
}
##################################Summary of NPE results
############################################
setwd(paste(root_path,"/NPE_FULL/NPE_Results",sep=""))
file_names <- dir()
d = length(file_names)
Record  = NULL
for(i in 1:d)
{
  selects = file_names[i]#BS 200 400
  split =  as.vector(str_split(selects,"_"))[[1]]
  M = split[3]
  covmatrix = split[5]
  n = as.numeric(split[7])
  p = as.numeric(split[9])
  model = as.numeric(split[11])
  delta = as.numeric(split[13])/10
  data = do.call(rbind,lapply(selects,read.csv)) 
  ERR = round(t(data)[2:8], 3)
  information = c(model,n,p,delta,M,covmatrix)
  ERR = c(information,ERR)
  Record = rbind(Record,ERR)
}
Record = data.frame(Record)
colnames(Record) = c("model","n","p","delta","M","covmatrix","Tn1","Size1","Size2","Size3","Time1","Time2","Time3")
#######################
setwd(paste(root_path,"/NPE_FULL/Summary",sep=""))
for(model in 1:5)
{
  for(covmatrix in c("AR","BCS"))
  {
    Record_model = NULL
    for(M in levels(factor(Record$M)))
    for(p in c(600,1200))
          {   
                Record_sub = Record[Record$model==model&Record$p==p&Record$covmatrix==covmatrix&Record$M==M,]
                Record_sub1 = Record_sub[order(as.numeric(Record_sub$delta)),]
                Record_model = rbind(Record_model,Record_sub1)
          }
    Record_model = as.data.frame(t(Record_model))
    write.csv(Record_model,paste("NPE_model_",model,"_covmatrix_",covmatrix,".csv",sep=""))
  }
}
###################################Summary of PPE results
############################################
setwd(paste(root_path,"/NPE_FULL/PPE_Results",sep=""))
file_names <- dir()
d = length(file_names)
Record  = NULL
for(i in 1:d)
{
  selects = file_names[i]#BS 200 400
  split =  as.vector(str_split(selects,"_"))[[1]]
  M = split[3]
  combine = split[5]
  n = as.numeric(split[7])
  p = as.numeric(split[9])
  model = as.numeric(split[11])
  delta = as.numeric(split[13])/10
  data = do.call(rbind,lapply(selects,read.csv)) 
  ERR = round(t(data)[2:5], 3)
  information = c(model,n,p,delta,M,combine)
  ERR = c(information,ERR)
  Record = rbind(Record,ERR)
}
Record = data.frame(Record)
colnames(Record) = c("model","n","p","delta","M","combine","Size1","Size2","Size3","Size_PPE")
#######################
setwd(paste(root_path,"/NPE_FULL/Summary",sep=""))
for(model in 5:5)
{
  for(combine in c("False","True"))
  {
    Record_model = NULL
      for(p in c(600))
      {   
        Record_sub = Record[Record$model==model&Record$p==p&Record$combine==combine,]
        Record_sub1 = Record_sub[order(as.numeric(Record_sub$delta)),]
        Record_model = rbind(Record_model,Record_sub1)
      }
    Record_model = as.data.frame(t(Record_model))
    write.csv(Record_model,paste("PPE_model_",model,"_combine_",combine,".csv",sep=""))
  }
}
