First Step : DATA Generation

---> Run 'DATA_GENERATE.m' to generate data.



Second Step : Fitting QUANTICO

---> Run 'AAA_QUANTICO_p_5_n_200.m', 'AAA_QUANTICO_p_10_n_200.m' and 'AAA_QUANTICO_p_20_n_200.m' for applying
     QUANTICO to (p=5,10,20;n=200) cases.



Third Step : Fitting VCQRM

---> Run 'AAA_VCQRM_p_5_n_200.m', 'AAA_VCQRM_p_10_n_200.m' and 'AAA_VCQRM_p_20_n_200.m'for applying VCQRM to
    (p=5,10,20;n=200) cases.



Fourth Step : Fitting LASSO-QR

---> Go to 'LASSO QR Method' folder.

---> Note that the whole generated dataset is copied inside a folder located in 'LASSO QR Method'folder, named 'Data'.

---> Run'QR_penalized_200_parallel.R' for applying LASSO-QR. User should change the the value of the variable
     'num_of_p' (set = 5,10,20) to obtain result for (p=5,10,20;n=200) cases.

---> Copy the generated csv files 'ZZZZ_MEAN_n_200_p_5_rep_25.csv', 'ZZZZ_MEAN_n_200_p_10_rep_25.csv', 'ZZZZ_MEAN_n_200_p_20_rep_25.csv',
     'ZZZZ_SD_n_200_p_5_rep_25.csv', 'ZZZZ_SD_n_200_p_10_rep_25.csv'  and 'ZZZZ_SD_n_200_p_20_rep_25.csv'to the parent folder 
    (i.e., 'Sample Size 200 Codes')


Fifth Step : Generating comparison table

---> Run 'Summarize.m' to combine all the results for 3 methods together.

---> The comparison tables are formed in the csv files (each dim : 7x9)
    (a) 'AAAA_WHOLE_SUMMARY_5_samplesize_200.csv'       (Case: Mean, p=5,n=200)
    (b) 'AAAA_WHOLE_SUMMARY_10_samplesize_200.csv'      (Case: Mean, p=10,n=200)
    (c) 'AAAA_WHOLE_SUMMARY_20_samplesize_200.csv'      (Case: Mean, p=20,n=200)
    (d) 'AAAA_WHOLE_SUMMARY_SD_5_samplesize_200.csv'    (Case: SD, p=5,n=200)
    (e) 'AAAA_WHOLE_SUMMARY_SD_10_samplesize_200.csv'   (Case: SD, p=10,n=200)
    (f) 'AAAA_WHOLE_SUMMARY_SD_20_samplesize_200.csv'   (Case: SD, p=20,n=200) 

---> The rows repersent : 
   (i)   pTPR
   (ii)  pFPR
   (iii) pAUC
   (iv)  gTPR
   (v)   gFPR
   (vi)  gAUC
   (vii) MSE

---> The columnss represent : 
   (i)   QUANTICO (\tau = 0.1)
   (ii)  VCQRM (\tau = 0.1)
   (iii) LASSO-QR (\tau = 0.1)
   (iv)  QUANTICO (\tau = 0.5)
   (v)   VCQRM (\tau = 0.5)
   (vi)  LASSO-QR (\tau = 0.5)
   (vii) QUANTICO (\tau = 0.9)
   (viii)VCQRM (\tau = 0.9)
   (ix)  LASSO-QR (\tau = 0.9)
   


---> Simulation results provided in the main paper is formed using these 4 csv files.
