DATA : 

---> 'SYNTHETIC_TCR.csv' contains synthetic data for 87 patients on 4 variables (Level 1 covariates) (NOT actual data)

---> 'SYNTHETIC_Mutation.csv' contains synthetic data for 87 patients on 7 variables (Level 2 covariates) (NOT actual data)

---> 'SYNTHETIC_CD8.csv' contains synthetic data for 87 patients on the Y variable (NOT actual data)

---> 'SYNTHETIC_clinical.csv' contains synthetic clinical data for 87 patients (NOT actual data)


QUANTICO ANALYSIS:

---> QUANTICO is performed on the synthetic dataset (imitating the real data) is
     performed in 'QUANTICO_SYNTHETIC_DATA_ANALYSIS.m'.

---> Note that, similar to what considered in the QUANTICO real data analysis,
     in the synthetic data also, since first 6 of the Level 2 covariates are often 
     zero and not varying much, in the QUANTICO analysis, we do NOT estimate the
     non-linear effect of those 6  Level 2 covariates. For them we only estimate
     the linear effect. However, for 7th Level 2 covariate, we estimate both 
     linear and non-linear effect.


Plots : 

---> Based on the output of the code 'QUANTICO_SYNTHETIC_DATA_ANALYSIS.m', running the code
     'QUANTICO_SYNTHETIC_DATA_ANALYSIS_PLOTS.R' will yield all the plots considered in the 
     main paper. Note that the results are diferent from the results in the main paper 
     since we use synthetic data here.

---> The patient-wise plot of mutation variables can be obtained using the code 
    'QUANTICO_SYNTHETIC_outlier_patient_mutations.R'.



