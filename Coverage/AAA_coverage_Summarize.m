

num_of_p = 10;
n = 100;

filename = ['ZZZ_All_coverage_',num2str(num_of_p),'_samplesize_',num2str(n),'_rep_50_iters_20000.csv'];
Coverage = csvread(filename);

filename2 = ['Unif_bound_width_p_',num2str(num_of_p),'_samplesize_',num2str(n),'_iters_20000.csv'];
width = csvread(filename2);

AAA= [width, 100*mean(Coverage(:,1)), 100*std(Coverage(:,1)),...
    1.5*sqrt(log(n))*width, 100*mean(Coverage(:,5)),100*std(Coverage(:,5))];

filename = ['AAA_coverage_SUMMARY_',num2str(num_of_p),'_',num2str(n),'.csv'];
csvwrite(filename,AAA) 

