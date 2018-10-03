function [Sig]=getcovcustom(sigma, lambda, lenrun, type)
if nargin<4
    type='gauss';
end
if strcmp(type, 'gauss')
    Sig=zeros(lenrun, lenrun);
    for ind=1:lenrun
        for jnd=1:lenrun
            Sig(ind, jnd)=sigma^2*exp(-.5*((ind-jnd)/lambda)^2);
        end
    end
elseif strcmp(type, 'ar1')
    Sig=zeros(lenrun, lenrun);
    for ind=1:lenrun
        for jnd=1:lenrun
            Sig(ind, jnd)=lambda^(abs(ind-jnd));
        end
    end
end
end