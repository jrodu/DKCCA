function [BigX, BigY, BigXnonoise, BigYnonoise, Acell, Bcell]=getsim(params)
% A sample params setup
lenrun=params.lenrun;
numtrials=params.numtrials;
numxelectrodes=params.numxelectrodes;
numyelectrodes=params.numyelectrodes;
startlagrange=params.startlagrange;
lag=params.lag;
ramplen=params.ramplen;
lenmax=params.lenmax;
totlenrho=lenmax+2*ramplen;
rho=params.rho;
arnoise=params.GaussnoiseCoeffs;
gaussnoise=params.GaussnoiseState;
gausslambdaCoeffs=params.GausslambdaCoeffs;
gausslambdaState=params.GausslambdaState;
%SignoiseA=getcovcustom(0, noiserho, numxelectrodes, 'ar1'); %  spatial
%SignoiseB=getcovcustom(0, noiserho, numyelectrodes, 'ar1'); % spatial
%SignoiseAR=getcovcustom(0, noiserho, lenrun, 'ar1')*.2;
% noise is gaussian process.  
SignoiseAR=getcovcustom(arnoise, gausslambdaCoeffs, lenrun); %not really AR here.
SignoiseGauss=getcovcustom(gaussnoise, gausslambdaState, lenrun);
%%
%parameters for coefficient matrix A

numparamsA=2;
Acell=cell(1, numparamsA);
muAcell=cell(1,numparamsA);
muAcell{1,1}=ones(1, lenrun);
muAcell{1,2}=zeros(1, lenrun);
sigmaA=1;
lambdaA=100;
for ind=1:numparamsA
    SigA=getcovcustom(sigmaA, lambdaA, lenrun);
    A= mvnrnd(muAcell{1, ind}, SigA, numxelectrodes);
    Acell{1, ind}=A;
end

%%
%parameters for coefficient matrix B

numparamsB=2;
Bcell=cell(1, numparamsB);
muBcell=cell(1,numparamsB);
muBcell{1,1}=ones(1, lenrun);
muBcell{1,2}=zeros(1, lenrun);
sigmaB=1;
lambdaB=100;
for ind=1:numparamsB
    SigB=getcovcustom(sigmaB, lambdaB, lenrun);
    B= mvnrnd(muBcell{1, ind}, SigB, numyelectrodes);
    Bcell{1, ind}=B;
end

%%
%B cov A
numparamsBcovA=1; %max is numparamsA
BcovAcell=cell(1, numparamsB);
muBcovAcell=cell(1,numparamsB);
muBcovAcell{1,1}=zeros(1, totlenrho);
sigmaBcovA=1;
lambdaBcovA=50;
for ind=1:numparamsBcovA
    SigB=getcovcustom(sigmaBcovA, lambdaBcovA, totlenrho);
    B= mvnrnd(muBcovAcell{1, ind}, SigB, numyelectrodes);
    BcovAcell{1, ind}=B;
end



%%
% covariance for latent state
lambdaX=params.lambdaX;
sigmaX=params.sigmaX;
SigX=getcovcustom(sigmaX, lambdaX, lenrun);

lambdaY=params.lambdaY;
sigmaY=params.sigmaY;
SigY=getcovcustom(sigmaY, lambdaY, lenrun);

%%
%simulate multiple trials
BigX=zeros(numxelectrodes, lenrun, numtrials);
BigXnonoise=zeros(size(BigX));
BigY=zeros(numyelectrodes, lenrun, numtrials);
BigYnonoise=zeros(size(BigY));
for ind=1:numtrials
    
    
    %%
    % simulate latent states for A
    xA=mvnrnd(zeros(1,lenrun), SigX, numparamsA);
    
    % and for B
    xB=mvnrnd(zeros(1,lenrun), SigY, numparamsB);
    
    %%
    % and simulate X
    X=Acell{1, 1}.*repmat(xA(1,:), numxelectrodes, 1);
    for jind=2:numparamsA
        X=X+Acell{1, jind}.*repmat(xA(jind,:), numxelectrodes, 1);
    end
    % and for Y
    Y=Bcell{1, 1}.*repmat(xB(1,:), numyelectrodes, 1);
    for jind=2:numparamsB
        Y=Y+Bcell{1, jind}.*repmat(xB(jind,:), numyelectrodes, 1);
    end
    
    Xnoise=mvnrnd(zeros(1, lenrun), SignoiseAR, numxelectrodes).*repmat(mvnrnd(zeros(1,lenrun), SignoiseGauss), numxelectrodes, 1);
    Ynoise=mvnrnd(zeros(1, lenrun), SignoiseAR, numyelectrodes).*repmat(mvnrnd(zeros(1,lenrun), SignoiseGauss), numyelectrodes, 1);

    
    %%
    startlag=datasample(startlagrange, 1);
    Astartlag=startlag-lag;
    
    YPrime=BcovAcell{1,1}.*repmat(ramprho(.8, 100, 20), numyelectrodes, 1).*repmat(ramprho(rho, lenmax, ramplen).*xA(1,Astartlag:Astartlag+totlenrho-1)+(1-ramprho(rho, lenmax, ramplen)).*mvnrnd(zeros(1,totlenrho), getcovcustom(sigmaBcovA, lambdaBcovA, totlenrho)), numyelectrodes, 1);
    %YPrime=BcovAcell{1,1}.*repmat(xA(1,Astartlag:Astartlag+totlenrho-1), numyelectrodes, 1)*.5;
    %comment out if you don't want correlation.
    if params.noiseflag
        Y(:, startlag:startlag+totlenrho-1)=Y(:, startlag:startlag+totlenrho-1).*repmat(1-ramprho(.8, 100, 20), numyelectrodes, 1)+YPrime;
    end
        %Y(:, startlag:startlag+totlenrho-1)=Y(:, startlag:startlag+totlenrho-1)*.5+YPrime;
    BigXnonoise(:,:,ind)=X;
    BigX(:,:,ind)=X+Xnoise;
    BigYnonoise(:,:,ind)=Y;
    BigY(:,:,ind)=Y+Ynoise;
end

end
