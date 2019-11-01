clear all

sample_size = 100;  % sample size, set = 100, 200 

num_of_p = 50;       % number of Level 1 covariates,
                     % generating 50 covariates, but used only 5,10,20
                     % in the studies

P = num_of_p+1;
G = 6;

QR_G_sim = ones(sample_size,G);
QR_P_sim = ones(sample_size,P);
QR_Y_sim = ones(sample_size,1);

rng(1)
taus = rand(sample_size,1);



for i = 1:sample_size
    QR_G_sim(i,:) = rand(1,G);
    QR_P_sim(i,2:P) = rand(1,P-1);
    QR_Y_sim(i) = data_y_1(taus(i),i,QR_P_sim,QR_G_sim);
end


filename1 = ['QR_G_sim1_',num2str(sample_size),'.csv'];
csvwrite(filename1,QR_G_sim)    
filename2 = ['QR_P_sim1_',num2str(sample_size),'.csv'];
csvwrite(filename2,QR_P_sim)  
filename3 = ['QR_Y_sim1_',num2str(sample_size),'.csv'];
csvwrite(filename3,QR_Y_sim)  

ts = linspace(0.05,0.95,19);

QR_Y_sim_TRUEs = ones(sample_size,length(ts));


for j = 1:length(ts)
    tau = ts(j);
    for i = 1:sample_size
        QR_Y_sim_TRUEs(i,j) = data_y_1(tau,i,QR_P_sim,QR_G_sim);
    end
end

csvwrite('QR_ts.csv',ts);
filename4 = ['QR_Y_sim1_TRUEs_',num2str(sample_size),'.csv'];
csvwrite(filename4,QR_Y_sim_TRUEs)  