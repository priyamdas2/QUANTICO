function ll=loglike(y,regressor,sigma,family)
    xb=sum(regressor,2);
    if family=='laplace' % (Laplace for quantile regression)
        nn = max(size(y));
        errors = y - xb;
        ll = nn*(log(sigma)+log(1-sigma));
        for ii = 1:nn
            if(errors(ii)>0)
                ll = ll - sigma*errors(ii);
            else
                ll = ll + (1-sigma)*errors(ii);
            end
        end
    end             
%     if family=='n'%normal w/ identity link
%         ll = -sum((y-xb).^2)/2/sigma-.5*log(sigma);
%     elseif family=='b'%binomial w/ logistic link
%         p = 1./(1+exp(-xb));
%         ll = sum(y.*log(p)+(1-y).*log(1-p));
%     elseif family=='p'%poisson w/ log link
%         lambda = exp(xb);
%         ll = sum(y.*log(lambda)-lambda-log(factorial(y)));
%     elseif family=='w'%weibull,w/ log link on scale parameter lambda, AFT, y(:,1) survival, y(:,2) censoring indicator
%         ll = sum(y(:,2).*(xb/sigma-log(sigma)+(1/sigma-1)*log(y(:,1))))-sum((y(:,1).*exp(xb)).^(1/sigma));
%     elseif family=='l'%log-normal w/ log link
%         ncens=sum(y(:,2));%number of non-censoring observations
%         n=length(y(:,2));
%         lfS = zeros(n,1);%log density and log survival
%         ind=(y(:,2)==1);
%         lfS(1:ncens)=-log(y(ind,1)*(sigma*2*pi)^.5)-(log(y(ind,1))-xb(ind)).^2/2/sigma;
%         lfS(ncens+1:end)=log(.5*erfc((log(y(~ind,1))-xb(~ind))/(2*sigma)^.5));
%         ll=sum(lfS);
%         %         ll=0;
%         %         for i = 1:length(xb)
%         %             ll = ll-lognlike([xb(i),sigma^.5],y(i,1),1-y(i,2));
%         %         end
%     end
end