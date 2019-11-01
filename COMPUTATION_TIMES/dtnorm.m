function d=dtnorm(x,mu,sig,a,b)
if nargin<5
    b = Inf;
    if nargin<4
        a = 0.01;
    end
end
alpha = (a-mu)/sig;
beta = (b-mu)/sig;

d = normpdf((x-mu)/sig)/sig/(normcdf(beta)-normcdf(alpha));