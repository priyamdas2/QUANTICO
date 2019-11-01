clear all

%%% G has dinesion 6, P has dimension 6 (we generate 51 but consider only
%%% first 6 while data analysis, among 6 first one is intercept)
sample_size = 200;
G = 6;
P = 50+1;

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

ts = linspace(0.1,0.9,9);

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