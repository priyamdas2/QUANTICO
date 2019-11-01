Frist Step : GENERATE DATA

---> Run 'AAA_DATA_GENERATE.m' to generate data for (n,p) = (100,5),(100,10),(200,5),(200,10),(200,20).




Second Step : ESTIMATE BAND WIDTH


---> Run 'AAA_QUANTICO_estimate_width_5_100.m' to estimate bandwith for the case (n,p) = (100,5).


---> Run 'AAA_QUANTICO_estimate_width_10_100.m' to estimate bandwith for the case (n,p) = (100,10).


---> Run 'AAA_QUANTICO_estimate_width_5_200.m' to estimate bandwith for the case (n,p) = (200,5).


---> Run 'AAA_QUANTICO_estimate_width_10_200.m' to estimate bandwith for the case (n,p) = (200,10).


---> Run 'AAA_QUANTICO_estimate_width_20_200.m' to estimate bandwith for the case (n,p) = (200,20).



Third Step : ESTIMATE COVERAGE (based on 50 replications)


---> Run 'AAA_QUANTICO_coverage_5_100.m' to estimate bandwith for the case (n,p) = (100,5).


---> Run 'AAA_QUANTICO_coverage_10_100.m' to estimate bandwith for the case (n,p) = (100,10).


---> Run 'AAA_QUANTICO_coverage_5_200.m' to estimate bandwith for the case (n,p) = (200,5).


---> Run 'AAA_QUANTICO_coverage_10_200.m' to estimate bandwith for the case (n,p) = (200,10).


---> Run 'AAA_QUANTICO_coverage_20_200.m' to estimate bandwith for the case (n,p) = (200,20). 



Fourth step : COMBINING RESULTS INTO TABLE


---> Run 'AAA_coverage_Summarize.m' for the settings (n,p) = (100,5),(100,10),(200,5),(200,10),(200,20).

---> Run 'AAA_coverage_FINAL_TABLE.m' to combine all case results into 1 table.

---> Generated table named 'AAA_coverage_FINAL_TABLE.csv' is of 5x6 dimension. This is reported in the main paper (or supplement).

---> The columns represent
   (i)   Bandwidth
   (ii)  Coverage (mean)
   (iii) Coverage (sd)
   (iv)  Inflated Bandwidth  (with inflation factor = 1.5*sqrt(n)
   (v)   Coverage for inflated band (mean)
   (vi)  Coverage for inflated band (sd)

---> The rows represent
   (i)   Case (n,p) = (100,5)
   (ii)  Case (n,p) = (100,10)
   (iii) Case (n,p) = (200,5)
   (iv)  Case (n,p) = (200,10)
   (v)   Case (n,p) = (200,20)





