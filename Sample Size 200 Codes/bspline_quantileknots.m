function B = bspline_quantileknots(n,nknots,x,y)
% B-spline basis function value matrix B(n) with quantile knots.
%
% Input arguments:
% n:
%    B-spline order (2 for linear, 3 for quadratic, etc.)
% nknots:
%    number of interior knots (total number of knots=nknots+2)
% x:
%    data points (could be matrix), used to decide the position of knots
% y(option):
%    points where the basis function is to be evaluated, has to have same
%    number of columns as x. If y is missing, y is set to x
%
% Output arguments:
% B:
%    Design matrix with dimension size(y,1) by size(y,2)*(n+nknots)
if nargin == 3
    y=x;
end
assert(size(x,2)==size(y,2),'x and y should have same number of columns');
G = size(x,2);
l = size(y,1);
M = nknots+n; %number of bases
B = zeros(l,G*M);
for g = 1:G
    x_tmp = x(:,g);
    xl = min(x_tmp); xr = max(x_tmp); %boundary knots
    knots.interior = quantileknots(x_tmp,nknots); %interior knots
    t = [xl*ones(n,1);knots.interior;xr*ones(n,1)]; %augmented knots
    B(:,(g-1)*M+1:g*M) = bspline_basismatrix(n,t,y(:,g));
end
