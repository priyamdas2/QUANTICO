
filename = ['AAA_coverage_SUMMARY_5_100.csv'];
file_1 = csvread(filename); 
filename = ['AAA_coverage_SUMMARY_10_100.csv'];
file_2 = csvread(filename); 
filename = ['AAA_coverage_SUMMARY_5_200.csv'];
file_3 = csvread(filename); 
filename = ['AAA_coverage_SUMMARY_10_200.csv'];
file_4 = csvread(filename); 
filename = ['AAA_coverage_SUMMARY_20_200.csv'];
file_5 = csvread(filename); 

FINAL = [file_1;file_2;file_3;file_4;file_5];

filename = ['AAA_coverage_FINAL_TABLE.csv'];
csvwrite(filename,FINAL);