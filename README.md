# High-Dim-Linear-Inference
 
Folder NPE_FULL

 -Functions.R This file contains all the necessary self-defined functions.
 
 -NPE_1.R The R script of the proposed method without power enhancement term JPE. The paper discusses numerous parameter settings; for simplicity in demonstration, we have set the parameters to fixed values.
 
 -PPE_1.R The R script of the partial penalized method in Shi et.al, (2019).
 
 -SIHR.R The R script of the SIHR package.

Due to the computational cost, we are providing all the outputs from the paper directly.

 NPE_Results: all the outputs of NPE.
 
 SIHR_Results: all the outputs of SIHR.
 
 PPE_Results: all the outputs of PPE

Sub-folder Summary: We use Data_Analysis.R to compile all the outputs into more easily readable Big Tables.

  -NPE_model_1_covmatrix_AR.csv: Model N1, covariance matrix Sigma is AR(1), the proposed method.
  
   delta: corresponding to different alternative hypothesis.
   
   Size1: empirical rejection rate of single sample splitting.
   
   Size2: empirical rejection rate of double sample splitting.
   
   Size3: empirical rejection rate of full sample.
   
   Time1, Time2, and Time3: the average computation time of single, double, and full sample method.

  -PPE_model_5_combine_False.csv: Model N5, single sample splitting, the partial penalized method.
  
   Size_PPE: empirical rejection rate of the PPE method.
   
   In single sample splitting, Size2 and Size3 are 0.
   
   In the double sample splitting "True", Size1 is 0.
   
  -SIHR_model_1_M_I.csv: The results of the SIHR method for model n1, using the SIHR-specified matrix

Folder PE_FULL

 -Functions.R This file contains all the necessary self-defined functions.
 
 -PE_1.R The R script of the proposed method with power enhancement term JPE. The paper discusses numerous parameter settings; for simplicity in demonstration, we have set the parameters to fixed values.
 
 Sub-folder Summary: We use Data_Analysis.R to compile all the outputs into more easily readable Big Tables.
 
  -ERR_model_1_combine_False.csv: The results of model 1 by using single sample splitting strategy.
  
   T_PE1: empirical rejection rate of single sample splitting strategy.
   
  -ERR_model_1_combine_True.csv: The results of model 1 by using double sample splitting strategy.
  
   T_PE2: empirical rejection rate of double sample splitting.
   
   T_PE3: empirical rejection rate of full sample.
   
   MCL: the empirical rejection rate of the MCL method. 
   
   Note for MCL method: Only in Models 1 and 2 under Method 1, there is a theoretical guarantee of the asymptotic property of the MCL method (Ma et al., 2021). Therefore, values from other models are excluded.

Folder Real Data Analysis

 -Real_Data_Analysis.R The R script to produce the output.
 
 -Real_functions.R All necessary functions.
 
 -HDGO.csv The IDs of selected group.
 
 -FULL_Y1_stat.csv: Outputs
 
  The first column is the GO ID.
  
  The second column is the value of proposed method without power enhancement term.
  
  The third column is the value of the method adding a power enhancement term.
  
  The fourth column is the value of ZC method.
  
 -transNOAH.Rdata: the real dataset.
   
Main Folder:

 Run Data_Analysis.R to summarize the results and save the output in the summary folder of each task folder.

