clear all

num_of_p = 5;
n = 100;

%%% Summary
filename = ['ZZZ_Mean_summary_',num2str(num_of_p),'_samplesize_',num2str(n),'_rep_24_iters_20000.csv'];
Mean_QUANTICO = csvread(filename);

filename = ['ZZZ_Sd_summary_',num2str(num_of_p),'_samplesize_',num2str(n),'_rep_24_iters_20000.csv'];
Sd_QUANTICO = csvread(filename);

filename = ['ZZZ_NO_CUT_Mean_summary_',num2str(num_of_p),'_samplesize_',num2str(n),'_rep_24_iters_20000.csv'];
Mean_VCQRM = csvread(filename);

filename = ['ZZZ_NO_CUT_Sd_summary_',num2str(num_of_p),'_samplesize_',num2str(n),'_rep_24_iters_20000.csv'];
Sd_VCQRM = csvread(filename);


filename = ['ZZZZ_MEAN_n_',num2str(n),'_p_',num2str(num_of_p),'_rep_25.csv'];
Mean_LASSOQR = csvread(filename);

filename = ['ZZZZ_SD_n_',num2str(n),'_p_',num2str(num_of_p),'_rep_25.csv'];
Sd_LASSOQR = csvread(filename);


WHOLE_SUMMARY = -.99*ones(7,9);

WHOLE_SUMMARY(:,[1,4,7]) = Mean_QUANTICO([1,2,7,3,4,8,10],[1,5,9]);
WHOLE_SUMMARY(:,[2,5,8]) = Mean_VCQRM([1,2,7,3,4,8,10],[1,5,9]);
WHOLE_SUMMARY([1,2,3,7],[3,6,9]) = Mean_LASSOQR(:,[1,5,9]);

WHOLE_SUMMARY = round(100*WHOLE_SUMMARY)/100;

WHOLE_SUMMARY_SD = -.99*ones(7,9);

WHOLE_SUMMARY_SD(:,[1,4,7]) = Sd_QUANTICO([1,2,7,3,4,8,10],[1,5,9]);
WHOLE_SUMMARY_SD(:,[2,5,8]) = Sd_VCQRM([1,2,7,3,4,8,10],[1,5,9]);
WHOLE_SUMMARY_SD([1,2,3,7],[3,6,9]) = Sd_LASSOQR(:,[1,5,9]);

WHOLE_SUMMARY_SD = round(100*WHOLE_SUMMARY_SD)/100;


filename_output = ['AAAA_WHOLE_SUMMARY_',num2str(num_of_p),'_samplesize_',num2str(n),'.csv'];
csvwrite(filename_output, WHOLE_SUMMARY);

filename_output = ['AAAA_WHOLE_SUMMARY_SD_',num2str(num_of_p),'_samplesize_',num2str(n),'.csv'];
csvwrite(filename_output, WHOLE_SUMMARY_SD);

