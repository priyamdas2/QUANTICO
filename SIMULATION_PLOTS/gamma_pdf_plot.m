


X=0:0.2:5; % Duration of data

Data=rand(1,100); % randomly generated 100 years of data

Z=gamfit(Data); % fitting distribution for shape and scale factor

alpha=4; % Shape factor
beta=0.25; % scale factor


X=0:0.2:5;
pdf=gampdf(X,alpha,beta); % pdf of the gamma distribution
plot(X,pdf,'-');