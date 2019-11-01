function knots = quantileknots(x,nknots)
% Get interior quantile knot sequence where the knots are the 
% quantiles from the empirical distribution of x
%
% Input:
% x: 
%    the data used to determine the quantile knots
% nknots:
%    number of interior knots
%
% Output:
% knots:
%    the interior quantile knots
x = unique(x);
n = length(x);
xsort = sort(x);
loc = n*(1:nknots)'./(nknots+1);
knots = xsort(round(loc));
