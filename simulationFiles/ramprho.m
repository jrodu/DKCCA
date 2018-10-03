function [rho]=ramprho(maxrho, lenmax, ramplen)
    rhoramp=(1:ramplen)*maxrho/ramplen;
    vecmaxrho=ones(1, lenmax)*maxrho;
    rho=[rhoramp vecmaxrho fliplr(rhoramp)];
end