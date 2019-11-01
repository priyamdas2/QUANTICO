clear all
tic;
ts = linspace(0.1,0.9,9);      % quantile levels
N = 50000;                     % number of MCMc iterations
thin = 5;
burnin = 10000;
post_size = (N-burnin)/thin;
family='laplace';


%%%%%%%%%%%%%%%%%%%% TCR (T) variable read %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P_raw = importdata('SYNTHETIC_TCR.csv');  % TCR_only
[n,p] = size(P_raw);
p = p+1;

P_var_names = [string('T.clonality');string('T.entropy');string('T.prod.cells');string('T.richness')];
%%%%%%%% PLOTS of TCR (T) variables (RAW)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
hold on
histogram(P_raw(:,1),'NumBins', 20)
xlabel('values', 'FontSize', 20);
ylabel('frequency', 'FontSize', 20);
title(['Synthetic T.clonality (raw)'])
hold off

figure
hold on
histogram(P_raw(:,2),'NumBins', 20)
xlabel('values', 'FontSize', 20);
ylabel('frequency', 'FontSize', 20);
title(['Synthetic T.entropy (raw)'])
hold off

figure
hold on
histogram(P_raw(:,3),'NumBins', 20)
xlabel('values', 'FontSize', 20);
ylabel('frequency', 'FontSize', 20);
title(['Synthetic T.prod.cells (raw)'])
hold off

figure
hold on
histogram(P_raw(:,4),'NumBins', 20)
xlabel('values', 'FontSize', 20);
ylabel('frequency', 'FontSize', 20);
title(['Synthetic T.richness (raw)'])
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%% Mutation (M) variable read %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G_raw = importdata('SYNTHETIC_Mutation.csv');
g = size(G_raw,2);

%%%%%%%%%%%%%%%%%%%% Immune (Y) variable read %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y_raw = importdata('SYNTHETIC_CD8.csv');


%%%%%%%%%%%%%%%%%%%% Immune (Y) transformation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y = log(y_raw);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%% TCR (P) transformation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P = ones(n,p);

mean_here = zeros(p,1);
sd_here = zeros(p,1);
for j = 2:p
    transformed = log(P_raw(:,j-1));
    mean_here(j) = sum(transformed)/max(size(transformed));
    sd_here(j) = std(transformed);
    P(:,j) = normcdf(transformed,mean_here(j),sd_here(j));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% PLOTS of TCR (T) variables (TRANSFORMED)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
hold on
histogram(P(:,2),'NumBins', 20)
xlabel('values', 'FontSize', 20);
ylabel('frequency', 'FontSize', 20);
title(['Synthetic T.clonality (transformed)'])
hold off

figure
hold on
histogram(P(:,3),'NumBins', 20)
xlabel('values', 'FontSize', 20);
ylabel('frequency', 'FontSize', 20);
title(['Synthetic T.entropy (transformed)'])
hold off

figure
hold on
histogram(P(:,4),'NumBins', 20)
xlabel('values', 'FontSize', 20);
ylabel('frequency', 'FontSize', 20);
title(['Synthetic T.prod.cells (transformed)'])
hold off

figure
hold on
histogram(P(:,5),'NumBins', 20)
xlabel('values', 'FontSize', 20);
ylabel('frequency', 'FontSize', 20);
title(['Synthetic T.richness (transformed)'])
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% Mutation (M) transformation (linear) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G = ones(n,g);
for k = 1:(g)
    G(:,k) = (G_raw(:,k) - (min(G_raw(:,k))))/((max(G_raw(:,k)))- (min(G_raw(:,k))));
end
rng(1)


%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% a_t and b_t controls cut off on parameter P
a_t = 4; b_t = 0.25;%gamma prior for t
% a_wt and b_w controls cut off on parameter P
a_w = (p-1)*(g-1)/2 +  (p-1);
b_w = (p-1)*(g-1)/2 +  (p-1);
a_tau = 5; b_tau = 25; 
s0 = 2.5*10^-4;
alpha_p = 0.1;
alpha_g = 0.1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

QR_Y_sim_ESTs = zeros(n,length(ts));
QR_Y_sim_ESTs_high = zeros(n,length(ts));
QR_Y_sim_ESTs_low = zeros(n,length(ts));

Prob_of_one = cell(length(ts),1); % for P
v_mean_all = cell(length(ts),1);
v_all = cell(length(ts),1);


Post_prob_linear = cell(length(ts),1);
Post_prob_nonlinear = cell(length(ts),1);
Post_prob_constanteffect = cell(length(ts),1);

G_selection_check = cell(length(ts),1);
aw_tracking = cell(length(ts),1);
bw_tracking = cell(length(ts),1);

P_post_indicators = zeros(n,p,length(ts));

betax_quantile_UP = cell(length(ts),p);
betax_quantile_DOWN = cell(length(ts),p);
G_post_indicators_linear = cell(length(ts),p);

betaX_quantile_UP = cell(length(ts),1);
betaX_quantile_DOWN = cell(length(ts),1);
G_post_indicators_nonlinear = cell(length(ts),1);
G_post_indicators_nonlinear_SIZED = cell(length(ts),1);



t_sd_vec = .1*[1/16,1/8,1/4,1/2,1,2,4,8,16,32,64,128,256,512,1024,2048];
eta_sd_vec = .1*[1/16,1/8,1/4,1/2,1,2,4,8,16,32,64,128,256,512,1024,2048];
xi_sd_vec=.1*[1/16,1/8,1/4,1/2,1,2,4,8,16,32,64,128,256,512,1024,2048];
summary_convergence = cell(length(ts),1);
patient_all_iter_estimates = cell(length(ts));

v_mean_array = zeros(n,p,length(ts));

v0_saver = cell(length(ts),1);
v_saver = cell(length(ts),1);
t_saver = cell(length(ts),1);


for jjj = 1:length(ts)
    rng(jjj)
    tau_now = jjj*0.1;
    t_old_array = zeros(N,1);
    jjj
    sigma = ts(jjj); % the quantile level
    [n,p]=size(P);
    x=cell(1,p);
    z=cell(1,p);
    for i = 1:p
        x{i}=G;
        z{i}=double.empty(n,0);
    end
    
    u = P;
    [n,p]=size(u);
    
    ximu_old = .1*randn(p,1);
    etamu_old = .1*randn(p,1);
    betamu = zeros(p,(N-burnin)/thin);
    betamu_old = etamu_old.*ximu_old;
    %design matrix
    inte = ones(n,1)/sqrt(n);%intercept
    %categorical
    Z = cell(p,1);
    Zind = cell(p,1);
    zmean = cell(p,1);
    zscale = cell(p,1);
    dZ = cell(p,1);
    xiZ_old = cell(p,1);
    betaZ = cell(p,1);
    betaZ_old = cell(p,1);
    etaZ_old = cell(p,1);
    nr=zeros(p,1);
    tauZ_old = cell(p,1);
    for j =1:p
        Z{j}=double.empty(n,0);
        Zind{j}=double.empty(0,1);
        zmean{j}=[];
        zscale{j}=[];
        nr(j) = size(z{j},2);
        for i =1:nr(j)
            try
                Ztmp = dummyvar(z{j}(:,i));
                Ztmp = Ztmp(:,2:end);
            catch
                Ztmp=z{j}(:,i);
            end
            zmean{j} = [zmean{j},mean(Ztmp)];
            Ztmp = Ztmp-repmat(mean(Ztmp),n,1);
            zscale{j} = [zscale{j},norm(Ztmp,'fro')*ones(1,size(Ztmp,2))];
            Ztmp = Ztmp/norm(Ztmp,'fro');
            Z{j} = [Z{j},Ztmp];
            Zind{j} = [Zind{j};i*ones(size(Ztmp,2),1)];
        end
        dZ{j}=size(Z{j},2);
        xiZ_old{j} = .1*randn(dZ{j},1);
        betaZ{j} = zeros(dZ{j},(N-burnin)/thin);
        etaZ_old{j}=.1*randn(nr(j),1);
        tauZ_old{j}=ones(nr(j),1);
        betaZ_old{j}=etaZ_old{j}(Zind{j}).*xiZ_old{j};
    end
    
    %continuous
    order = 4; %order of splines
    nknots = 16; %number of interior knots
    M = order + nknots;
    K = makeK(M);
    Kinv = pinv(K);
    dX = cell(p,1);
    X=cell(p,1);
    Xind=cell(p,1);
    xix_old = cell(p,1);
    xiX_old=cell(p,1);
    betaX = cell(p,1);
    betaX_old=cell(p,1);
    aa_X = cell(p,1);
    etax_old = cell(p,1);
    etaX_old = cell(p,1);
    xx= cell(p,1);
    %xxt=cell(p,1);
    betax = cell(p,1);
    betax_old = cell(p,1);
    aa_x = cell(p,1);
    nq = zeros(p,1);
    taux_old = cell(p,1);
    tauX_old = cell(p,1);
    for j = 1:p
        % Only g-th G variable is continuos, so for that variable only we
        % use B-spline. Since the rest of the coordinates are binary, we
        % don't use B-spline.
        xx{j}=bspline_quantileknots(order,nknots,x{j}(:,g),x{j}(:,g));
        %     if ifpredu==1
        %         xxt{j}=bspline_quantileknots(order,nknots,x{j},xt{j});
        %     end
        X{j} = double.empty(n,0);
        Xind{j} =double.empty(0,1);
        nq(j) = size(x{j},2); % (= g in this case)
        for i = nq(j):nq(j)
            Xtmp = xx{j}; % xx{j}(:,(i-1)*M+1:i*M);
            [U,S,~]=svd(Xtmp*Kinv*Xtmp','econ');
            S = diag(S);
            nullvals = S < 10^-10;
            d = max(3, find( cumsum(S(~nullvals))/sum(S(~nullvals)) > .995 ,1));
            d = min(d, sum(~nullvals));
            Xtmp = U(:,1:d)*diag(sqrt(S(1:d)));
            Xtmp2 = [ones(n,1),x{j}(:,i)];
            Xtmp = Xtmp- Xtmp2*(Xtmp2\Xtmp);
            Xtmp = Xtmp/norm(Xtmp,'fro');
            X{j} = Xtmp; %[X{j},Xtmp];
            % Xind{j} = [Xind{j};i*ones(size(Xtmp,2),1)];
        end
        Xind{j} = ones(size(Xtmp,2),1);
        
        for i = 1:nq(j)
            x{j}(:,i) = x{j}(:,i)-mean(x{j}(:,i));
            x{j}(:,i) = x{j}(:,i)/norm(x{j}(:,i),'fro');
        end
        
        dX{j}= size(X{j},2);
        xix_old{j} = .1*randn(nq(j),1);
        xiX_old{j} = .1*randn(dX{j},1);
        betaX{j} = zeros(dX{j},(N-burnin)/thin);
        etax_old{j} = .1*randn(nq(j),1);
        etaX_old{j} = .1*randn(1,1); % since only 1 cont G %.1*randn(nq(j),1);
        betax{j}=zeros(nq(j),(N-burnin)/thin);
        betax_old{j} = etax_old{j}.*xix_old{j};
        taux_old{j}=ones(nq(j),1);
        tauX_old{j}=ones(1,1);
        betaX_old{j} = etaX_old{j}(Xind{j}).*xiX_old{j};
    end
    
    
    
    
    for j = 1:p
        for i = 1:nq(j)
            xibar=mean(abs(xix_old{j}(i)),1);
            etax_old{j}(i) = etax_old{j}(i).*xibar;
            xix_old{j}(i) = xix_old{j}(i)./xibar;
            if(i == nq(j))
                xibar=mean(abs(xiX_old{j}),1);
                etaX_old{j} = etaX_old{j}.*xibar;
                xiX_old{j} = xiX_old{j}./repmat(xibar,[sum(Xind{j}==1),1]);
            end
        end
        for i = 1:nr(j)
            xibar=mean(abs(xiZ_old{j}(Zind{j}==i)),1);
            etaZ_old{j}(i) = etaZ_old{j}(i).*xibar;
            xiZ_old{j}(Zind{j}==i) = xiZ_old{j}(Zind{j}==i)./repmat(xibar,[sum(Zind{j}==i),1]);
            
        end
    end
    etamu_old = etamu_old.*abs(ximu_old);
    ximu_old = ximu_old./abs(ximu_old);
    t = zeros((N-burnin)/thin,1);
    ll = zeros((N-burnin)/thin,1);
    t_old=t(1);
    
    
    gammax=cell(p,1);
    gammaX=cell(p,1);
    gammaZ=cell(p,1);
    gammax_old=cell(p,1);
    gammaX_old=cell(p,1);
    gammaZ_old=cell(p,1);
    
    gammax{1} = ones(nq(1),(N-burnin)/thin);
    gammaX{1} = ones(1,(N-burnin)/thin); %ones(nq(1),(N-burnin)/thin);
    gammaZ{1} = ones(nr(1),(N-burnin)/thin);
    gammax_old{1} = gammax{1}(:,1);
    gammaX_old{1} = gammaX{1}(:,1);
    gammaZ_old{1} = gammaZ{1}(:,1);
    
    for j =2:p
        gammax{j} = s0*ones(nq(j),(N-burnin)/thin);
        gammaX{j} =  s0*ones(1,(N-burnin)/thin); % s0*ones(nq(j),(N-burnin)/thin);
        gammaZ{j} = s0*ones(nr(j),(N-burnin)/thin);
        gammax_old{j} = gammax{j}(:,1);
        gammaX_old{j} = gammaX{j}(:,1);
        gammaZ_old{j} = gammaZ{j}(:,1);
    end
    gammamu = s0*ones(p,(N-burnin)/thin);
    gammamu_old = gammamu(:,1);
    v = zeros(n,p,(N-burnin)/thin);
    v0 = zeros(n,p,(N-burnin)/thin);
    v_saver{jjj} = zeros(n,p,(N-burnin)/thin);
    v0_saver{jjj} = zeros(n,p,(N-burnin)/thin);
    t_saver{jjj} = zeros((N-burnin)/thin,1);
    
    ac_eta=0;
    ac_xi =0;
    ac_t=0;
    eta_iter=5;
    xi_iter=5;
    t_iter=5;
    eta_sd=eta_sd_vec(eta_iter);
    xi_sd=xi_sd_vec(xi_iter);
    t_sd=t_sd_vec(t_iter);
    sigma_old = sigma;
    taumu_old = ones(p,1);
    omega_old = .2;
    regressor_old=zeros(n,p);
    for j = 1:p
        regressor_old(:,j) = u(:,j).*threshold(betamu_old(j)*inte+x{j}*betax_old{j}+X{j}*betaX_old{j}+Z{j}*betaZ_old{j},t_old);
    end
    %algorithm begins
    %tic;
    ll_old = loglike(y,regressor_old,sigma_old,family);
    iter=0;
    for mc = 2:N
        %update eta
        
        for j = 1:p
            for i = 1:nq(j)
                regressor_new = regressor_old;
                etax_new = etax_old{j};
                etax_new(i) = etax_old{j}(i) + eta_sd*randn(1);
                betax_new = etax_new.*xix_old{j};
                regressor_new(:,j) = u(:,j).*threshold(betamu_old(j)*inte+x{j}*betax_new+X{j}*betaX_old{j}+Z{j}*betaZ_old{j},t_old);
                ll_new = loglike(y,regressor_new,sigma_old,family);
                lr = ll_new-ll_old-.5*(etax_new(i).^2./gammax_old{j}(i)./taux_old{j}(i))...
                    +.5*(etax_old{j}(i).^2./gammax_old{j}(i)./taux_old{j}(i));
                if lr > log(rand(1))
                    ac_eta = ac_eta+1;
                    ll_old=ll_new;
                    etax_old{j}(i) = etax_new(i);
                    betax_old{j} = betax_new;
                    regressor_old(:,j) = regressor_new(:,j);
                end
            end
        end
        for j = 1:p
            for i = nq(j):nq(j)
                regressor_new = regressor_old;
                etaX_new = etaX_old{j} + eta_sd*randn(1);
                betaX_new = etaX_new.*xiX_old{j};
                regressor_new(:,j) = u(:,j).*threshold(betamu_old(j)*inte+x{j}*betax_old{j}+X{j}*betaX_new+Z{j}*betaZ_old{j},t_old);
                ll_new = loglike(y,regressor_new,sigma_old,family);
                lr = ll_new-ll_old-.5*(etaX_new.^2./gammaX_old{j}./tauX_old{j})...
                    +.5*(etaX_old{j}.^2./gammaX_old{j}./tauX_old{j});
                if lr > log(rand(1))
                    ac_eta = ac_eta+1;
                    ll_old=ll_new;
                    etaX_old{j} = etaX_new;
                    betaX_old{j} = betaX_new;
                    regressor_old(:,j) = regressor_new(:,j);
                end
            end
        end
        for j = 1:p
            for i = 1:nr(j)
                regressor_new = regressor_old;
                etaZ_new = etaZ_old{j};
                etaZ_new(i) = etaZ_old{j}(i) + eta_sd*randn(1);
                betaZ_new = etaZ_new(Zind{j}).*xiZ_old{j};
                regressor_new(:,j) = u(:,j).*threshold(betamu_old(j)*inte+x{j}*betax_old{j}+X{j}*betaX_old{j}+Z{j}*betaZ_new,t_old);
                ll_new = loglike(y,regressor_new,sigma_old,family);
                lr = ll_new-ll_old-.5*(etaZ_new(i).^2./gammaZ_old{j}(i)./tauZ_old{j}(i))...
                    +.5*(etaZ_old{j}(i).^2./gammaZ_old{j}(i)./tauZ_old{j}(i));
                if lr > log(rand(1))
                    ac_eta = ac_eta+1;
                    ll_old=ll_new;
                    etaZ_old{j}(i) = etaZ_new(i);
                    betaZ_old{j} = betaZ_new;
                    regressor_old(:,j) = regressor_new(:,j);
                end
            end
        end
        for j = 1:p
            regressor_new = regressor_old;
            etamu_new = etamu_old(j) + eta_sd*randn(1);
            betamu_new = etamu_new.*ximu_old(j);
            regressor_new(:,j) = u(:,j).*threshold(betamu_new*inte+x{j}*betax_old{j}+X{j}*betaX_old{j}+Z{j}*betaZ_old{j},t_old);
            ll_new = loglike(y,regressor_new,sigma_old,family);
            lr = ll_new-ll_old-.5*(etamu_new.^2./gammamu_old(j)./taumu_old(j))...
                +.5*(etamu_old(j).^2./gammamu_old(j)./taumu_old(j));
            if lr > log(rand(1))
                ac_eta = ac_eta+1;
                ll_old=ll_new;
                etamu_old(j) = etamu_new;
                betamu_old(j) = betamu_new;
                regressor_old(:,j) = regressor_new(:,j);
            end
        end
        
        %update mmu,mx,mX,mZ CHANGED on 01/08/2019
        
        mmu_old=(2*binornd(1,1./(1+exp(-2*ximu_old)))-1);
        mx_old = cell(p,1);
        mX_old = cell(p,1);
        mZ_old = cell(p,1);
        
        %update xi
        for j = 1:p
            mx_old{j}=(2*binornd(1,1./(1+exp(-2*xix_old{j})))-1);
            mX_old{j}=(2*binornd(1,1./(1+exp(-2*xiX_old{j})))-1);
            mZ_old{j}=(2*binornd(1,1./(1+exp(-2*xiZ_old{j})))-1);
            regressor_new = regressor_old;
            xix_new = xix_old{j} + xi_sd*randn(nq(j),1);
            betax_new = etax_old{j}.*xix_new;
            xiX_new = xiX_old{j} + xi_sd*randn(dX{j},1);
            betaX_new = etaX_old{j}(Xind{j}).*xiX_new;
            xiZ_new = xiZ_old{j} + xi_sd*randn(dZ{j},1);
            betaZ_new = etaZ_old{j}(Zind{j}).*xiZ_new;
            ximu_new = ximu_old(j) + xi_sd*randn(1);
            betamu_new = etamu_old(j).*ximu_new;
            regressor_new(:,j) = u(:,j).*threshold(betamu_new*inte+x{j}*betax_new+X{j}*betaX_new+Z{j}*betaZ_new,t_old);
            ll_new = loglike(y,regressor_new,sigma_old,family);
            lr = ll_new-ll_old-.5*(sum(sum((xix_new-mx_old{j}).^2))+sum(sum((xiX_new-mX_old{j}).^2))+sum(sum((xiZ_new-mZ_old{j}).^2))+sum((ximu_new-mmu_old(j)).^2))...
                +.5*(sum(sum((xix_old{j}-mx_old{j}).^2))+sum(sum((xiX_old{j}-mX_old{j}).^2))+sum(sum((xiZ_old{j}-mZ_old{j}).^2))+sum((ximu_old(j)-mmu_old(j)).^2));
            if lr > log(rand(1))
                ac_xi = ac_xi+1;
                ll_old=ll_new;
                xix_old{j} = xix_new;
                betax_old{j} = betax_new;
                xiX_old{j} = xiX_new;
                betaX_old{j} = betaX_new;
                xiZ_old{j} = xiZ_new;
                betaZ_old{j} = betaZ_new;
                ximu_old(j) = ximu_new;
                betamu_old(j) = betamu_new;
                regressor_old(:,j) = regressor_new(:,j);
            end
        end
        
        
        %rescale eta and xi
        for j =1:p
            for i = 1:nq(j)
                xibar=abs(xix_old{j}(i));
                etax_old{j}(i) = etax_old{j}(i).*xibar;
                xix_old{j}(i) = xix_old{j}(i)./xibar;
                if(i == nq(j))
                    xibar=mean(abs(xiX_old{j}),1);
                    etaX_old{j} = etaX_old{j}.*xibar;
                    xiX_old{j} = xiX_old{j}./repmat(xibar,[sum(Xind{j}==1),1]);
                end
            end
            for i = 1:nr(j)
                xibar=mean(abs(xiZ_old{j}(Zind{j}==i)),1);
                etaZ_old{j}(i) = etaZ_old{j}(i).*xibar;
                xiZ_old{j}(Zind{j}==i) = xiZ_old{j}(Zind{j}==i)./repmat(xibar,[sum(Zind{j}==i),1]);
            end
        end
        etamu_old = etamu_old.*abs(ximu_old);
        ximu_old = ximu_old./abs(ximu_old);
        
        
        %update t
        t_old_array(mc) = t_old;
        t_new = t_old + normrnd(0,t_sd);
        for j = 1:p
            regressor_new(:,j) = u(:,j).*threshold(betamu_old(j)*inte+x{j}*betax_old{j}+X{j}*betaX_old{j}+Z{j}*betaZ_old{j},t_new);
        end
        ll_new = loglike(y,regressor_new,sigma_old,family);
        lr = ll_new-ll_old + log(gampdf(t_new,a_t,b_t))- log(gampdf(t_old,a_t,b_t));
        if lr>log(rand(1))
            ac_t = ac_t + 1;
            t_old = t_new;
            regressor_old = regressor_new;
            ll_old=ll_new;
        end
        
        %update tau
        for j =1:p
            taux_old{j}=1./gamrnd(a_tau+.5,1./(b_tau+etax_old{j}.^2./gammax_old{j}/2));
            tauX_old{j}=1./gamrnd(a_tau+.5,1./(b_tau+etaX_old{j}.^2./gammaX_old{j}/2));
            tauZ_old{j}=1./gamrnd(a_tau+.5,1./(b_tau+etaZ_old{j}.^2./gammaZ_old{j}/2));
        end
        taumu_old=1./gamrnd(a_tau+.5,1./(b_tau+etamu_old.^2./gammamu_old/2));
        
        
        %update gammma
        aw = a_w;
        bw = b_w;
        for j = 2:p  % new edition
            ptmp=sqrt(s0)*omega_old/(1-omega_old)*exp((1-s0)*etax_old{j}.^2/2/s0./taux_old{j});
            gammax_old{j}=binornd(1,ptmp./(1+ptmp));
            gammax_old{j}(ptmp==Inf)=1;
            gammax_old{j}(gammax_old{j}==0)=s0;
            ptmp=sqrt(s0)*omega_old/(1-omega_old)*exp((1-s0)*etaX_old{j}.^2/2/s0./tauX_old{j});
            gammaX_old{j}=binornd(1,ptmp./(1+ptmp));
            gammaX_old{j}(ptmp==Inf)=1;
            gammaX_old{j}(gammaX_old{j}==0)=s0;
            ptmp=sqrt(s0)*omega_old/(1-omega_old)*exp((1-s0)*etaZ_old{j}.^2/2/s0./tauZ_old{j});
            gammaZ_old{j}=binornd(1,ptmp./(1+ptmp));
            gammaZ_old{j}(ptmp==Inf)=1;
            gammaZ_old{j}(gammaZ_old{j}==0)=s0;
            ngamma=sum(gammax_old{j}==1)+sum(gammaX_old{j}==1)+sum(gammaZ_old{j}==1);
            aw = aw+ngamma;
            bw = bw+2*nq(j)-(nq(j)-1)+nr(j)-ngamma; % (since (nq(j)-1) only has linear part)
        end
        ptmp=sqrt(s0)*omega_old/(1-omega_old)*exp((1-s0)*etamu_old.^2/2/s0./taumu_old);
        gammamu_old=binornd(1,ptmp./(1+ptmp));
        gammamu_old(ptmp==Inf)=1;
        gammamu_old(gammamu_old==0)=s0;
        gammamu_old(1) = 1;  % new edition
        ngamma = sum(gammamu_old(2:end)==1);
        aw = aw+ngamma;
        bw = bw+p-ngamma-1;
        %update omega (rho)
        omega_old = betarnd(aw,bw);
        aw_tracking{jjj}(mc) = aw;
        bw_tracking{jjj}(mc) = bw;
        
        
        
        if mod(mc,1000)==0
            fprintf('%d th quantile...%d steps finished, %d steps to go\n', tau_now, mc, N-mc)
        end
        
        if mod(mc,100)==0&&mc<=burnin
            ac_eta_mornitor = ac_eta/(mc-1)/(2*sum(nq)+sum(nr)+p);
            ac_xi_mornitor = ac_xi/(mc-1)/p;
            ac_t_mornitor = ac_t/(mc-1);
            if ac_eta_mornitor>.5&&eta_iter<9
                eta_iter=eta_iter+1;
                eta_sd=eta_sd_vec(eta_iter);
            elseif ac_eta_mornitor<.2&&eta_iter>1
                eta_iter=eta_iter-1;
                eta_sd=eta_sd_vec(eta_iter);
            end
            if ac_xi_mornitor>.5&&xi_iter<9
                xi_iter=xi_iter+1;
                xi_sd=xi_sd_vec(xi_iter);
            elseif ac_xi_mornitor<.2&&xi_iter>1
                xi_iter=xi_iter-1;
                xi_sd=xi_sd_vec(xi_iter);
            end
            if ac_t_mornitor>.5&&t_iter<9
                t_iter=t_iter+1;
                t_sd=t_sd_vec(t_iter);
            elseif ac_t_mornitor<.2&&t_iter>1
                t_iter=t_iter-1;
                t_sd=t_sd_vec(t_iter);
            end
        end
        
        if(mc == burnin)
            summary_convergence{jjj} = [ac_eta_mornitor,ac_xi_mornitor,ac_t_mornitor;eta_sd,xi_sd,t_sd];
        end
        
        %save samples every "thin" iterations
        if mod(mc,thin)==0&&mc>burnin
            iter = iter + 1;
            for j =1:p
                betax{j}(:,iter)=betax_old{j};
                betaX{j}(:,iter)=betaX_old{j};
                betaZ{j}(:,iter)=betaZ_old{j};
                v0(:,j,iter) = betamu_old(j)*(gammamu_old(j)==1)*inte+...
                    x{j}*(betax_old{j}.*(gammax_old{j}==1))+...
                    X{j}*(betaX_old{j}.*(gammaX_old{j}(Xind{j})==1))+...
                    Z{j}*(betaZ_old{j}.*(gammaZ_old{j}(Zind{j})==1));
                
                aa_x{j}(:,iter) = (betax_old{j}.*(gammax_old{j}==1));
                aa_X{j}(:,iter) = (betaX_old{j}.*(gammaX_old{j}(Xind{j})==1));
                gammax{j}(:,iter)=gammax_old{j};
                gammaX{j}(:,iter)=gammaX_old{j};
                gammaZ{j}(:,iter)=gammaZ_old{j};
            end
            betamu(:,iter)=betamu_old;
            t(iter) = t_old;
            ll(iter) = ll_old;
            gammamu(:,iter)=gammamu_old;
            %v(:,:,iter) = threshold(v0(:,:,iter),t_old);
            v(:,:,iter) = threshold(v0(:,:,iter),t(iter));
            v(:,1,iter) = v0(:,1,iter);
            v0_saver{jjj}(:,:,iter) = v0(:,:,iter);
            v_saver{jjj}(:,:,iter) = v(:,:,iter);
            t_saver{jjj}(iter) = t(iter);
        end
    end
    
    %%% Posterior inference P
    Prob_of_one{jjj} = zeros(n,p);
    v_mean = mean(v(:,:,1:post_size),3);
    v_mean_all{jjj} = v_mean;
    v_all{jjj} = v;
    
    v_quantile_UP = zeros(n,p);
    v_quantile_DOWN = zeros(n,p);
    
    for ii = 1:n
        for jj = 1:p
            v_quantile_UP(ii,jj) = quantile(v(ii,jj,1:post_size),1-alpha_p/2);
            v_quantile_DOWN(ii,jj) = quantile(v(ii,jj,1:post_size),alpha_p/2);
        end
    end
    P_post_indicators(:,:,jjj) = (sign(v_quantile_UP).*sign(v_quantile_DOWN)==1);
    
    patient_all_iter_estimates{jjj} = zeros(n,post_size);
    for iter = 1:post_size
        patient_all_iter_estimates{jjj}(:,iter) = sum(P.*v(:,:,iter),2);
    end
    
    for i = 1:n
        for j = 1:p
            count = 0;
            for iter = 1:post_size
                if(abs(v(i,j,iter)) < 10^(-10))
                    count = count + 1;
                end
            end
            Prob_of_one{jjj}(i,j) = (post_size - count)/post_size;
            if(Prob_of_one{jjj}(i,j) < 0.5)
                v_mean(i,j) = 0;
            end
        end
    end
    
    %%% Posterior inference G
    Post_prob_linear{jjj} = zeros(g,p);
    Post_prob_nonlinear{jjj} = zeros(1,p);
    Post_prob_constanteffect{jjj} = zeros(p,1);
    
    
    for j =1:p
        Post_prob_linear{jjj}(:,j) = sum((gammax{j}==1),2)/iter;
        Post_prob_nonlinear{jjj}(:,j) = sum((gammaX{j}==1))/iter;
        Post_prob_constanteffect{jjj}(j) = sum((gammamu(j,:)==1),2)/iter;
    end
    %G_selection_check{jjj} = max((Post_prob_linear{jjj}>0.5),(Post_prob_nonlinear{jjj}>0.5));
    
    
    %     % linear
    %     betax_quantile_UP{jjj} = zeros(nq(j),p);
    %     betax_quantile_DOWN{jjj} = zeros(nq(j),p);
    %     for j = 1:p
    %         for k = 1:nq(j)
    %             betax_quantile_UP{jjj}(k,j) = quantile(aa_x{j}(k,1:post_size),1-alpha_g/2);
    %             betax_quantile_DOWN{jjj}(k,j) = quantile(aa_x{j}(k,1:post_size),alpha_g/2);
    %         end
    %         G_post_indicators_linear{jjj}(:,j) = (abs(sign(betax_quantile_UP{jjj}(:,j))+sign(betax_quantile_DOWN{jjj}(:,j)))>=1);
    %     end
    %
    %     % non-linear
    %     betaX_quantile_UP{jjj} = zeros(dX{j},p);
    %     betaX_quantile_DOWN{jjj} = zeros(dX{j},p);
    %     for j = 1:p
    %         for k = 1:dX{j}
    %             betaX_quantile_UP{jjj}(k,j) = quantile(aa_X{j}(k,1:post_size),1-alpha_g/2);
    %             betaX_quantile_DOWN{jjj}(k,j) = quantile(aa_X{j}(k,1:post_size),alpha_g/2);
    %         end
    %         G_post_indicators_nonlinear{jjj}(:,j) = (abs(sign(betaX_quantile_UP{jjj}(:,j))+sign(betaX_quantile_DOWN{jjj}(:,j)))>=1);
    %         for k = 1:nq(j)
    %             G_post_indicators_nonlinear_SIZED{jjj}(k,j) = max(G_post_indicators_nonlinear{jjj}((Xind{j}==k),j));
    %         end
    %     end
    
    
    betax{j}(:,iter);
    
    sum(abs(v_mean)>0.00001)
    QR_Y_sim_ESTs_high(:,jjj) = transpose(quantile(transpose(patient_all_iter_estimates{jjj}),0.975));
    QR_Y_sim_ESTs(:,jjj) = sum(P.*v_mean,2);
    QR_Y_sim_ESTs_low(:,jjj) = transpose(quantile(transpose(patient_all_iter_estimates{jjj}),0.025));
    v_mean_array(:,:,jjj) = v_mean;
    
    
    
end


%%% G selection using FDR cutoff %%%%%%

num_of_cuts = 201;
cutoffs = linspace(0,1,num_of_cuts);
cutoffs = [0.05,0.1,0.2,0.8, 0.9,0.95];

Post_prob = cell(length(ts),1);
Post_prob_arrays = cell(length(ts),1);
Post_prob_q_value_arrays = cell(length(ts),1);
aa = cell(length(ts),1);
bb = cell(length(ts),1);
cumsum_by_position = cell(length(ts),1);
selected_G = cell(length(ts),length(cutoffs));





for jjj = 1:length(ts)
    Post_prob{jjj} = zeros(g,p-1);
    Post_prob_arrays{jjj} = zeros(g*(p-1),1);
    
    for jj = 1:(p-1)
        for ii = 1:g
            if(ii == g)
                Post_prob{jjj}(ii,jj) = max(Post_prob_nonlinear{jjj}(1,jj+1), Post_prob_linear{jjj}(ii,jj+1));
            else
                Post_prob{jjj}(ii,jj) = Post_prob_linear{jjj}(ii,jj+1);
            end
            Post_prob_arrays{jjj}((ii-1)*(p-1)+jj) = Post_prob{jjj}(ii,jj);
        end
    end
    Post_prob_q_value_arrays{jjj} = 1 - Post_prob_arrays{jjj};
    [aa{jjj},bb{jjj}] = sort(Post_prob_q_value_arrays{jjj});
    
    cum_sum_aa = cumsum(aa{jjj});
    cumsum_by_position{jjj} = zeros(g*(p-1),1);
    for i = 1:(g*(p-1))
        cumsum_by_position{jjj}(i) = cum_sum_aa(i)/i;
    end
    for cut = 1:length(cutoffs)
        %%%% G variable TP/FP %%%%
        positions_selected = bb{jjj}(cumsum_by_position{jjj}<cutoffs(cut));
        %Post_prob_q_value_arrays{jjj}(positions_to_be_included)
        selected_G_array = zeros(g*(p-1),1);
        selected_G_array(positions_selected) = 1;
        
        selected_G{jjj,cut} = zeros(g,p-1);
        for ii = 1:g
            for jj = 1:(p-1)
                selected_G{jjj,cut}(ii,jj)= selected_G_array((ii-1)*(p-1)+jj);
            end
        end
        
    end
    
end

% for i = 1:10
%
%     figure
%     hold on
%     h = zeros(2,1);
%     h(1) = plot(ts, QR_Y_sim_ESTs(i,:),'color','b');
%     h(2) = plot(ts, QR_Y_sim_ESTs_high(i,:),'--','color','g');
%     h(3) = plot(ts, QR_Y_sim_ESTs_low(i,:),'--','color','g');
%     axis([0 1 -3 3])
%     legend([h(1)],{'BEHAVIOR_{QR}'},'Location','northwest');
%     title(['Estimated Quantile values for Patient number = ' num2str(i)])
%     xlabel('\tau')
%     ylabel('Y')
%     hold off
%
% end


%%% P posterior inference
patient_nums = 1:n;
Mean_probs = zeros(length(ts),p);
Mean_probs_patient_T_prod_cells = zeros(length(ts), n);

for jjj = 1:length(ts)
    Mean_probs(jjj,:) = mean(Prob_of_one{jjj},1);
    Mean_probs_patient_T_prod_cells(jjj,:) = transpose(Prob_of_one{jjj}(:,4));
end

%%% G posterior inference


Mean_probs_G = cell(p,1);
for j = 1:p
    Mean_probs_G{j} = zeros(size(G,2),length(ts));
    for jjj = 1:length(ts)
        Mean_probs_G{j}(1:(size(G,2)-1),jjj) = Post_prob_linear{jjj}(1:(size(G,2)-1),j);
        Mean_probs_G{j}(size(G,2),jjj) = max(Post_prob_nonlinear{jjj}(1,j),Post_prob_linear{jjj}(size(G,2),j));
    end
end

Mean_probs_G_constant = zeros(p,length(ts));

for jjj = 1:length(ts)
    for j = 1:p
        Mean_probs_G_constant(j,jjj) = Post_prob_constanteffect{jjj}(j);
    end
end





% figure
% hold on
% col=hsv(p);
% h = zeros(p-1,1);
% for i = 2:p
%     h(i-1) = plot(ts(1,:),transpose(Mean_probs(:,i)),'color',col(i,:));
% end
% axis([0 1 0 1])
% legend([h(1) h(2) h(3) h(4)],...
%     {'T.clonality', 'T.entropy', ...
%     'T.prod.cells', 'T.richness'},'Location','southwest');
% line([0,1], [0.5, 0.5], 'Color', 'g','LineStyle',':');
% title(['Posterior probabilities of having non-zero coefficients'])
% xlabel('Quantiles')
% ylabel('Posterior probabilities')
% hold off







%%% G selection plots


% which_P = 1;
%
% figure
% hold on
% col=hsv(g);
% h = zeros(g,1);
% for i = 1:g
%     h(i) = plot(ts(1,:),transpose(Mean_probs_G{which_P+1}(i,:)),'color',col(i,:));
% end
% axis([0 1 0 1])
% legend([h(7) h(1) h(2) h(3) h(4) h(5) h(6)],{'TMB','TTN','MUC16', 'RYR2', 'CSMD3', 'TP53', 'USH2A'},'Location','northwest');
% line([0,1], [0.5, 0.5], 'Color', 'g','LineStyle',':');
% title(['Posterior probabilities of having non-zero coefficients G for ',P_var_names(which_P)'])
% xlabel('Quantiles')
% ylabel('Posterior probabilities')
% hold off
%


% which_P = 2;
%
% figure
% hold on
% col=hsv(g);
% h = zeros(g,1);
% for i = 1:g
%     h(i) = plot(ts(1,:),transpose(Mean_probs_G{which_P+1}(i,:)),'color',col(i,:));
% end
% axis([0 1 0 1])
% legend([h(7) h(1) h(2) h(3) h(4) h(5) h(6)],{'TMB','TTN','MUC16', 'RYR2', 'CSMD3', 'TP53', 'USH2A'},'Location','northwest');
% line([0,1], [0.5, 0.5], 'Color', 'g','LineStyle',':');
% title(['Posterior probabilities of having non-zero coefficients G for ',P_var_names(which_P)'])
% xlabel('Quantiles')
% ylabel('Posterior probabilities')
% hold off


% which_P = 3;
% 
% figure
% hold on
% col=hsv(g);
% h = zeros(g,1);
% for i = 1:g
%     h(i) = plot(ts(1,:),transpose(Mean_probs_G{which_P+1}(i,:)),'color',col(i,:));
% end
% axis([0 1 0 1])
% legend([h(7) h(1) h(2) h(3) h(4) h(5) h(6)],{'TMB','TTN','MUC16', 'RYR2', 'CSMD3', 'TP53', 'USH2A'},'Location','northwest');
% line([0,1], [0.5, 0.5], 'Color', 'g','LineStyle',':');
% title(['Posterior probabilities of having non-zero coefficients G for ',P_var_names(which_P)'])
% xlabel('Quantiles')
% ylabel('Posterior probabilities')
% hold off



% which_P = 4;
%
% figure
% hold on
% col=hsv(g);
% h = zeros(g,1);
% for i = 1:g
%     h(i) = plot(ts(1,:),transpose(Mean_probs_G{which_P+1}(i,:)),'color',col(i,:));
% end
% axis([0 1 0 1])
% legend([h(7) h(1) h(2) h(3) h(4) h(5) h(6)],{'TMB','TTN','MUC16', 'RYR2', 'CSMD3', 'TP53', 'USH2A'},'Location','northwest');
% line([0,1], [0.5, 0.5], 'Color', 'g','LineStyle',':');
% title(['Posterior probabilities of having non-zero coefficients G for ',P_var_names(which_P)'])
% xlabel('Quantiles')
% ylabel('Posterior probabilities')
% hold off
%

%%%%%%%%%%%%%% beta_plot %%%%%%%%%%%%%%%


mean_coef = zeros(length(ts), p);
mean_coef_patient_wise_prod_cells = zeros(length(ts), n);
mean_coef_UP = zeros(length(ts), p);
mean_coef_DOWN = zeros(length(ts), p);

for jjj = 1:length(ts)
    mean_coef(jjj,:) = mean(v_mean_all{jjj},1);
    mean_coef_patient_wise_prod_cells(jjj,:) = transpose(v_mean_all{jjj}(:,4));
    mean_now = zeros(p,post_size);
    for ii = 1:p
        for iter = 1:post_size
            elements = zeros(n,1);
            for kk = 1:n
                elements(kk) = v_all{jjj}(kk,ii,iter);
            end
            mean_now(ii,iter) = mean(elements);
        end
    end
    for ii = 1:p
        mean_coef_UP(jjj,ii) = quantile(mean_now(ii,:), 0.975);
        mean_coef_DOWN(jjj,ii) = quantile(mean_now(ii,:), 0.025);
    end
end

% which_P_beta = 4;

% figure
% hold on
% col=hsv(3);
% h(1) = plot(ts(1,:),transpose(mean_coef(:,which_P_beta)/sd_here(which_P_beta)),'color','b');
% h(2) = plot(ts(1,:),transpose(mean_coef_UP(:,which_P_beta)/sd_here(which_P_beta)),'--','color','r');
% h(3) = plot(ts(1,:),transpose(mean_coef_DOWN(:,which_P_beta)/sd_here(which_P_beta)),'--','color','r');
% axis([0 1 -4 12])
% title(['Posterior mean of coefficient of ',P_var_names(which_P)'])
% xlabel('Quantiles')
% ylabel('Posterior mean of coefficient of T.prod.cells')
% hold off


filename_output = ['AAA_to_R_mean_probs_patient_T_prod_cells_iters_',num2str(N),'.csv'];
csvwrite(filename_output, Mean_probs_patient_T_prod_cells);

filename_output = ['AAA_to_R_mean_coef_patient_wise_prod_cells_iters_',num2str(N),'.csv'];
csvwrite(filename_output, mean_coef_patient_wise_prod_cells);



filename_output = ['AAA_to_R_Mean_probs_G_TCR_prod_cells_iters_',num2str(N),'.csv'];
csvwrite(filename_output, Mean_probs_G{3+1});

filename_output = ['AAA_to_R_Mean_probs_P_iters_',num2str(N),'.csv'];
csvwrite(filename_output, Mean_probs);

