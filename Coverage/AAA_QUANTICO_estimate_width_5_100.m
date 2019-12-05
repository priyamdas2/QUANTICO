clear all
num_of_p = 5;
sample_size = 100;
ts = linspace(0.05,0.95,19);
N = 20000;
thin = 5;
burnin = 5000;
post_size = (N-burnin)/thin;
%cutoff=0.5;
family='laplace';
filename1 = ['QR_G_sim1_',num2str(sample_size),'.csv'];    
filename2 = ['QR_P_sim1_',num2str(sample_size),'.csv'];
filename3 = ['QR_Y_sim1_',num2str(sample_size),'.csv'];
filename4 = ['QR_Y_sim1_TRUEs_',num2str(sample_size),'.csv'];
G = csvread(filename1);
P_all = csvread(filename2);
P = P_all(:,1:(1+num_of_p));
y_raw = csvread(filename3);
FDR_eta = 0.2;

QR_Y_sim_TRUEs_raw = csvread(filename4);

y = y_raw;

QR_Y_sim_TRUEs = QR_Y_sim_TRUEs_raw;


[n,p]=size(P);
g = size(G,2);
%parameter
% a_t and b_t controls cut off on parameter P
a_t = 4; b_t = 0.25;%gamma prior for t


% a_wt and b_w controls cut off on parameter P
a_w = num_of_p*g;b_w=num_of_p*g;

a_tau = 5; b_tau = 25;
s0 = 2.5*10^-4;
%a_w = 1/2;b_w=1/2;
%a_w = 2*5+2*5+5;b_w= 2*5+2*5+5;

alpha_p = 0.1;
alpha_g = 0.1;


QR_Y_sim_ESTs = zeros(n,length(ts));
QR_Y_sim_ESTs_high = zeros(n,length(ts));
QR_Y_sim_ESTs_low = zeros(n,length(ts));

Prob_of_one = cell(length(ts),1); % for p
v_mean_all = cell(length(ts),1);


Post_prob_linear = cell(length(ts),1);
Post_prob_nonlinear = cell(length(ts),1);

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



t_sd_vec = .1*[1/16,1/8,1/4,1/2,1,2,4,8,16];
eta_sd_vec = .1*[1/16,1/8,1/4,1/2,1,2,4,8,16];
xi_sd_vec=.1*[1/16,1/8,1/4,1/2,1,2,4,8,16];
summary_convergence = cell(length(ts),1);
patient_all_iter_estimates = cell(length(ts));
parfor jjj = 1:length(ts)
    rng(jjj+50)
    tau_now = 0.05*jjj;
    t_old_array = zeros(N,1);
    sigma = ts(jjj); % the quantile level
    [n,p]=size(P);
    x=cell(1,p);
    z=cell(1,p);
    for i = 1:p
        x{i}=G;
        z{i}=double.empty(n,0);
    end 
    
    
%     rel = version('-release');
%     if (all(rel>='2011a'))
%         rng('default');
%         rng(seed);
%     end
    
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
        xx{j}=bspline_quantileknots(order,nknots,x{j},x{j});
        %     if ifpredu==1
        %         xxt{j}=bspline_quantileknots(order,nknots,x{j},xt{j});
        %     end
        X{j} = double.empty(n,0);
        Xind{j} =double.empty(0,1);
        nq(j) = size(x{j},2);
        for i = 1:nq(j)
            Xtmp = xx{j}(:,(i-1)*M+1:i*M);
            [U,S,~]=svd(Xtmp*Kinv*Xtmp','econ');
            S = diag(S);
            nullvals = S < 10^-10;
            d = max(3, find( cumsum(S(~nullvals))/sum(S(~nullvals)) > .995 ,1));
            d = min(d, sum(~nullvals));
            Xtmp = U(:,1:d)*diag(sqrt(S(1:d)));
            Xtmp2 = [ones(n,1),x{j}(:,i)];
            Xtmp = Xtmp- Xtmp2*(Xtmp2\Xtmp);
            Xtmp = Xtmp/norm(Xtmp,'fro');
            x{j}(:,i) = x{j}(:,i)-mean(x{j}(:,i));
            x{j}(:,i) = x{j}(:,i)/norm(x{j}(:,i),'fro');
            X{j} = [X{j},Xtmp];
            Xind{j} = [Xind{j};i*ones(size(Xtmp,2),1)];
        end
        dX{j}=size(X{j},2);
        xix_old{j} = .1*randn(nq(j),1);
        xiX_old{j} = .1*randn(dX{j},1);
        betaX{j} = zeros(dX{j},(N-burnin)/thin);
        etax_old{j} = .1*randn(nq(j),1);
        etaX_old{j} = .1*randn(nq(j),1);
        betax{j}=zeros(nq(j),(N-burnin)/thin);
        betax_old{j} = etax_old{j}.*xix_old{j};
        taux_old{j}=ones(nq(j),1);
        tauX_old{j}=ones(nq(j),1);
        betaX_old{j} = etaX_old{j}(Xind{j}).*xiX_old{j};
    end
    
    
    
    
    for j = 1:p
        for i = 1:nq(j)
            xibar=mean(abs(xix_old{j}(i)),1);
            etax_old{j}(i) = etax_old{j}(i).*xibar;
            xix_old{j}(i) = xix_old{j}(i)./xibar;
            xibar=mean(abs(xiX_old{j}(Xind{j}==i)),1);
            etaX_old{j}(i) = etaX_old{j}(i).*xibar;
            xiX_old{j}(Xind{j}==i) = xiX_old{j}(Xind{j}==i)./repmat(xibar,[sum(Xind{j}==i),1]);
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
    gammaX{1} = ones(nq(1),(N-burnin)/thin);
    gammaZ{1} = ones(nr(1),(N-burnin)/thin);
    gammax_old{1} = gammax{1}(:,1);
    gammaX_old{1} = gammaX{1}(:,1);
    gammaZ_old{1} = gammaZ{1}(:,1);
    
    for j =2:p
        gammax{j} = s0*ones(nq(j),(N-burnin)/thin);
        gammaX{j} = s0*ones(nq(j),(N-burnin)/thin);
        gammaZ{j} = s0*ones(nr(j),(N-burnin)/thin);
        gammax_old{j} = gammax{j}(:,1);
        gammaX_old{j} = gammaX{j}(:,1);
        gammaZ_old{j} = gammaZ{j}(:,1);
    end
    gammamu = s0*ones(p,(N-burnin)/thin);
    gammamu_old = gammamu(:,1);
    v = zeros(n,p,(N-burnin)/thin); 
    v0 = zeros(n,p,(N-burnin)/thin);
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
            for i = 1:nq(j)
                regressor_new = regressor_old;
                etaX_new = etaX_old{j};
                etaX_new(i) = etaX_old{j}(i) + eta_sd*randn(1);
                betaX_new = etaX_new(Xind{j}).*xiX_old{j};
                regressor_new(:,j) = u(:,j).*threshold(betamu_old(j)*inte+x{j}*betax_old{j}+X{j}*betaX_new+Z{j}*betaZ_old{j},t_old);
                ll_new = loglike(y,regressor_new,sigma_old,family);
                lr = ll_new-ll_old-.5*(etaX_new(i).^2./gammaX_old{j}(i)./tauX_old{j}(i))...
                    +.5*(etaX_old{j}(i).^2./gammaX_old{j}(i)./tauX_old{j}(i));
                if lr > log(rand(1))
                    ac_eta = ac_eta+1;
                    ll_old=ll_new;
                    etaX_old{j}(i) = etaX_new(i);
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
        
        %update
        
        mmu_old=2*binornd(1,1./(1+exp(-2*ximu_old)))-1;
        mx_old = cell(p,1);
        mX_old = cell(p,1);
        mZ_old = cell(p,1);
        
        %update xi
        for j = 1:p
            mx_old{j}=2*binornd(1,1./(1+exp(-2*xix_old{j})))-1;
            mX_old{j}=2*binornd(1,1./(1+exp(-2*xiX_old{j})))-1;
            mZ_old{j}=2*binornd(1,1./(1+exp(-2*xiZ_old{j})))-1;
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
                xibar=mean(abs(xiX_old{j}(Xind{j}==i)),1);
                etaX_old{j}(i) = etaX_old{j}(i).*xibar;
                xiX_old{j}(Xind{j}==i) = xiX_old{j}(Xind{j}==i)./repmat(xibar,[sum(Xind{j}==i),1]);
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
            bw = bw+2*nq(j)+nr(j)-ngamma;
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
        end
    end
    

    %%% Posterior inference P
    Prob_of_one{jjj} = zeros(n,p);
    v_mean = mean(v(:,:,1:post_size),3);
    v_mean_all{jjj} = v_mean;
    
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
    for ff = 1:post_size
        patient_all_iter_estimates{jjj}(:,ff) = sum(P.*v(:,:,ff),2);
    end
    
    for i = 1:n
        for j = 1:p
            count = 0;
            for ff = 1:post_size
                if(abs(v(i,j,ff)) < 10^(-10))
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
    Post_prob_nonlinear{jjj} = zeros(g,p);
    
    
    for j =1:p
        Post_prob_linear{jjj}(:,j) = sum((gammax{j}==1),2)/iter;
        Post_prob_nonlinear{jjj}(:,j) = sum((gammaX{j}==1),2)/iter;
    end
    %G_selection_check{jjj} = max((Post_prob_linear{jjj}>0.5),(Post_prob_nonlinear{jjj}>0.5));
    
    
    % linear
    betax_quantile_UP{jjj} = zeros(nq(j),p);
    betax_quantile_DOWN{jjj} = zeros(nq(j),p);
    for j = 1:p
        for k = 1:nq(j)
            betax_quantile_UP{jjj}(k,j) = quantile(aa_x{j}(k,1:post_size),1-alpha_g/2);
            betax_quantile_DOWN{jjj}(k,j) = quantile(aa_x{j}(k,1:post_size),alpha_g/2);
        end
        G_post_indicators_linear{jjj}(:,j) = (abs(sign(betax_quantile_UP{jjj}(:,j))+sign(betax_quantile_DOWN{jjj}(:,j)))>=1);
    end
    
    % non-linear
    betaX_quantile_UP{jjj} = zeros(dX{j},p);
    betaX_quantile_DOWN{jjj} = zeros(dX{j},p);
    for j = 1:p
        for k = 1:dX{j}
            betaX_quantile_UP{jjj}(k,j) = quantile(aa_X{j}(k,1:post_size),1-alpha_g/2);
            betaX_quantile_DOWN{jjj}(k,j) = quantile(aa_X{j}(k,1:post_size),alpha_g/2);
        end
        G_post_indicators_nonlinear{jjj}(:,j) = (abs(sign(betaX_quantile_UP{jjj}(:,j))+sign(betaX_quantile_DOWN{jjj}(:,j)))>=1);
        for k = 1:nq(j)
            G_post_indicators_nonlinear_SIZED{jjj}(k,j) = max(G_post_indicators_nonlinear{jjj}((Xind{j}==k),j));
        end
    end
    
    
    betax{j}(:,iter)
    
    sum(abs(v_mean)>0.00001)
    QR_Y_sim_ESTs_high(:,jjj) = transpose(quantile(transpose(patient_all_iter_estimates{jjj}),0.975));
    QR_Y_sim_ESTs(:,jjj) = sum(P.*v_mean,2);
    QR_Y_sim_ESTs_low(:,jjj) = transpose(quantile(transpose(patient_all_iter_estimates{jjj}),0.025));
end


%%% AUC calculation %%%%%%%%%%%%%%



null_coords_P = cell(length(ts),1);
true_coords_P = cell(length(ts),1);
for jjj = 1:5
    null_coords_P{jjj} = 4:p;
    true_coords_P{jjj} = [2,3];
end

for jjj = 6:length(ts)
    null_coords_P{jjj} = 5:p;
    true_coords_P{jjj} = [2,3,4];
end

total_non_zero_G = 2;
Post_prob = cell(length(ts),1);
Post_prob_arrays = cell(length(ts),1);
cumsum_by_position = cell(length(ts),1);
Post_prob_q_value_arrays = cell(length(ts),1);
aa = cell(length(ts),1);
bb = cell(length(ts),1);


num_of_cuts = 201;
cutoffs = linspace(0,1,num_of_cuts);
TP_P_AUC = zeros(length(ts),length(cutoffs));
FP_P_AUC = zeros(length(ts),length(cutoffs));

TP_G_FDR_AUC = zeros(length(ts),length(cutoffs));
FP_G_FDR_AUC = zeros(length(ts),length(cutoffs));

TPR_total = cell(length(ts),1); 
FPR_total = cell(length(ts),1);
TPR_FPR_total = cell(length(ts),1);

AUC = zeros(length(ts),1);



for jjj = 1:length(ts)
    Post_prob{jjj} = zeros(g,p-1);
    for ii = 1:g
        for jj = 1:(p-1)
            Post_prob{jjj}(ii,jj) = max(Post_prob_nonlinear{jjj}(ii,jj+1), Post_prob_linear{jjj}(ii,jj+1));
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
        
        %%%% P variable TP/FP %%%%
        TP_P_AUC(jjj,cut) = sum(sum(Prob_of_one{jjj}(:,true_coords_P{jjj})>cutoffs(cut)));
        FP_P_AUC(jjj,cut) = sum(sum(Prob_of_one{jjj}(:,null_coords_P{jjj})>cutoffs(cut)));
        
        %%%% G variable TP/FP %%%%
        positions_selected = bb{jjj}(cumsum_by_position{jjj}<cutoffs(cut));
        %Post_prob_q_value_arrays{jjj}(positions_to_be_included)
        selected_G_array_AUC = zeros(g*(p-1),1);
        selected_G_array_AUC(positions_selected) = 1;
        
        selected_G_AUC = zeros(g,p-1);
        for ii = 1:g
            for jj = 1:(p-1)
                selected_G_AUC(ii,jj)= selected_G_array_AUC((ii-1)*(p-1)+jj);
            end
        end
        
        TP_G_FDR_AUC(jjj,cut) = selected_G_AUC(1,2-1) + selected_G_AUC(2,3-1);
        FP_G_FDR_AUC(jjj,cut) = sum(sum(selected_G_AUC)) - TP_G_FDR_AUC(jjj,cut);
    end
    
    TPR_total{jjj}= zeros(length(cutoffs),length(cutoffs));
    FPR_total{jjj}= zeros(length(cutoffs),length(cutoffs));
    TPR_FPR_total{jjj}= zeros(length(cutoffs)^2,2);
    
    for i_cut = 1:length(cutoffs)
        for j_cut = 1:length(cutoffs)
            TPR_total{jjj}(i_cut,j_cut) = (TP_P_AUC(jjj,i_cut) + TP_G_FDR_AUC(jjj,j_cut))/...
                (n*max(size(true_coords_P{jjj}))+ total_non_zero_G);
            FPR_total{jjj}(i_cut,j_cut) = (FP_P_AUC(jjj,i_cut) + FP_G_FDR_AUC(jjj,j_cut))/...
                (n*max(size(null_coords_P{jjj}))+ g*(p-1) - total_non_zero_G);
            TPR_FPR_total{jjj}((i_cut-1)*length(cutoffs)+j_cut,1) = TPR_total{jjj}(i_cut,j_cut);
            TPR_FPR_total{jjj}((i_cut-1)*length(cutoffs)+j_cut,2) = FPR_total{jjj}(i_cut,j_cut);
        end
    end
    uniques_TPR_FPR_total = unique(TPR_FPR_total{jjj},'rows');
    
    [aaa,bbb] = sort(uniques_TPR_FPR_total(:,2));
    sorted_mat = zeros(max(size(uniques_TPR_FPR_total)),2);
    sorted_mat(:,2) = aaa;
    sorted_mat(:,1) = uniques_TPR_FPR_total(bbb,1);
    AUC(jjj) = trapz(sorted_mat(:,2),sorted_mat(:,1));
end


%%% TRUE POSITIVE RATE (P) %%%

TPR_P = zeros(length(ts),1);
FPR_P = zeros(length(ts),1);

for jjj = 1:length(ts)
    TPR_P(jjj) = (1/(max(size(true_coords_P{jjj}))*n))*sum(sum(Prob_of_one{jjj}(:,true_coords_P{jjj})>0.5));%.*P_post_indicators(:,true_coords_P,jjj)>0.5));
    FPR_P(jjj) = (1/(max(size(null_coords_P{jjj}))*n))*sum(sum(Prob_of_one{jjj}(:,null_coords_P{jjj})>0.5));%.*P_post_indicators(:,null_coords_P,jjj)>0.5));
end


%%% TRUE POSITIVE RATE (G) %%%

total_non_zero_G = 2;
Post_prob = cell(length(ts),1);
Post_prob_arrays = cell(length(ts),1);
selected_G_matrix = cell(length(ts),1);
Post_prob_q_value_arrays = cell(length(ts),1);
TPR_FPR_FDR = zeros(2,length(ts));
for jjj = 1:length(ts)
        %%% FDR   
    positions_selected = bb{jjj}(cumsum_by_position{jjj}<=FDR_eta);
    %Post_prob_q_value_arrays{jjj}(positions_to_be_included)
    selected_G_array = zeros(g*(p-1),1);
    selected_G_array(positions_selected) = 1;
    
    selected_G_matrix{jjj} = zeros(g,p-1);
    for ii = 1:g
        for jj = 1:(p-1)
            selected_G_matrix{jjj}(ii,jj)= selected_G_array((ii-1)*(p-1)+jj);
        end
    end
    
    TP_all_FDR = selected_G_matrix{jjj}(1,2-1) + selected_G_matrix{jjj}(2,3-1);
    FP_all_FDR = sum(sum(selected_G_matrix{jjj})) - TP_all_FDR;
    TPR_FPR_FDR(1,jjj) = TP_all_FDR/total_non_zero_G;
    TPR_FPR_FDR(2,jjj) = FP_all_FDR/(g*(p-1) - total_non_zero_G);
end

%%% MCC and Combined TPR-FPR

Mean_probs_G = cell(p,1);
for jjj = 1:length(ts)
    for j = 1:p
        Mean_probs_G{j}(:,jjj) = max(Post_prob_nonlinear{jjj}(:,j),Post_prob_linear{jjj}(:,j));
    end
end


Total_TP = zeros(length(ts),1);
Total_FP = zeros(length(ts),1);
Total_TN = zeros(length(ts),1);
Total_FN = zeros(length(ts),1);
Total_TPR_FPR = zeros(2,length(ts));
MCC = zeros(length(ts),1);

for jjj = 1:length(ts)
    Total_TP(jjj) = TPR_P(jjj)*max(size(true_coords_P{jjj}))*n + TPR_FPR_FDR(1,jjj)*total_non_zero_G;
    Total_FP(jjj) = FPR_P(jjj)*(p-1-max(size(true_coords_P{jjj})))*n + TPR_FPR_FDR(2,jjj)*(g*(p-1) - total_non_zero_G);
    Total_TN(jjj) = max(size(true_coords_P{jjj}))*n + total_non_zero_G - Total_TP(jjj);
    Total_FN(jjj) = (p-1-max(size(true_coords_P{jjj})))*n + (g*(p-1) - total_non_zero_G) - Total_FP(jjj);
    
    Total_TPR_FPR(1,jjj) = Total_TP(jjj)/(max(size(true_coords_P{jjj}))*n + total_non_zero_G);
    Total_TPR_FPR(2,jjj) = Total_FP(jjj)/((p-1-max(size(true_coords_P{jjj})))*n + (g*(p-1) - total_non_zero_G));
    
    MCC(jjj) = (Total_TP(jjj)*Total_TN(jjj) - Total_FP(jjj)*Total_FN(jjj))/sqrt((Total_TP(jjj) + Total_FP(jjj))*...
        (Total_TP(jjj)+Total_FN(jjj))*(Total_TN(jjj) + Total_FP(jjj))*(Total_TN(jjj) + Total_FN(jjj)));
end

%%% Patient specific uniform bound

patient_estimate_diff_max = abs(patient_all_iter_estimates{1} - QR_Y_sim_ESTs(:,1));

for jjj = 2:length(ts)
    patient_estimate_diff_max = max(patient_estimate_diff_max,abs(patient_all_iter_estimates{jjj}-QR_Y_sim_ESTs(:,jjj)));
end

patient_estimate_diff_max_array = reshape(patient_estimate_diff_max,... 
[size(patient_estimate_diff_max,1)*size(patient_estimate_diff_max,2),1]);


unif_bound_width = quantile(patient_estimate_diff_max_array,0.95);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




filename_output_1 = ['Unif_bound_width_p_',num2str(num_of_p),'_samplesize_',num2str(sample_size),'_iters_',num2str(N),'.csv']; 
csvwrite(filename_output_1,unif_bound_width);