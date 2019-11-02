clear all
tic;
sample_size = 200;
n = sample_size;
g = 6;

no_simulations = 25;


ts = linspace(0.1,0.9,9);
for zzzz = 1:no_simulations
    rng(zzzz+100)
    %%%% Data Generate
    p = 50+1;
    QR_G_sim = ones(sample_size,g);
    QR_P_sim = ones(sample_size,p);
    QR_Y_sim = ones(sample_size,1);
    
    taus = rand(sample_size,1);
    for i = 1:sample_size
        QR_G_sim(i,:) = rand(1,g);
        QR_P_sim(i,2:p) = rand(1,p-1);
        QR_Y_sim(i) = data_y_1(taus(i),i,QR_P_sim,QR_G_sim);
    end
    
    QR_Y_sim_TRUEs = zeros(sample_size, length(ts));
    for j = 1:length(ts)
        tau = ts(j);
        for i = 1:sample_size
            QR_Y_sim_TRUEs(i,j) = data_y_1(tau,i,QR_P_sim,QR_G_sim);
        end
    end
    
    filename_output = ['QR_Y_sim_TRUEs_n_',num2str(sample_size),'_seed_',num2str(zzzz),'.csv'];
    csvwrite(filename_output, QR_Y_sim_TRUEs);
    
    filename_output = ['QR_G_sim_n_',num2str(sample_size),'_seed_',num2str(zzzz),'.csv'];
    csvwrite(filename_output, QR_G_sim);
    
    filename_output = ['QR_P_sim_n_',num2str(sample_size),'_seed_',num2str(zzzz),'.csv'];
    csvwrite(filename_output, QR_P_sim);
    
    filename_output = ['QR_Y_sim_n_',num2str(sample_size),'_seed_',num2str(zzzz),'.csv'];
    csvwrite(filename_output, QR_Y_sim);
    
end