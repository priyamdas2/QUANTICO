function x=rtnorm(mu,sig,a,b)
if nargin<4
    b = Inf;
    if nargin<3
        a = 0.01;
    end
end
alpha = (a-mu)/sig;
beta = (b-mu)/sig;
cdfalpha = normcdf(alpha);
x = norminv(cdfalpha+rand(1)*(normcdf(beta)-cdfalpha),0,1)*sig+mu;
    