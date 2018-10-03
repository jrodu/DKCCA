addpath('.\simulationFiles', '.\util');

sig=1; %toggle this parameter to set the noise level in the simulations
%       this value correspnds to params.GaussnoiseState below

%settings for DKCCA algorithm
numpermute=500;  % this is too low, but is set here for speed.  
setnumperwindow=1;
windowlength=21;
numcomponents=1;
regwindow=.7:.01:1;
gsdval=5;


%parameters for the simulation
params.lenrun=500;
params.numtrials=180;
params.numxelectrodes=96;
params.numyelectrodes=16;
params.startlagrange=290:300;
params.GaussnoiseCoeffs=1;
params.GaussnoiseState=sig;
params.lag=20;
params.ramplen=20;
params.lenmax=100;
params.rho=.8;
params.sigmaX=1;
params.sigmaY=1;
params.lambdaX=40;
params.lambdaY=20;
params.GausslambdaCoeffs=30;
params.GausslambdaState=80;
params.noiseflag=1;

%simulate data
datacell=cell(2, 1);
[datacell{1,1}, datacell{2,1}]=getsim(params);


%run DKCCA algorithm
[dcell]=MDCCG(datacell{1,1}, datacell{2,1}, setnumperwindow, windowlength, numpermute, numcomponents, regwindow, gsdval);


%image plots

subplot(1,2,1)
imagesc(dcell.original)
caxis([0 1])
axis square
titleuse=sprintf("pre-inference \\sigma_{noise}=%d", sig);
title(titleuse, 'Interpreter', 'tex')
xlabel('time (ms), region Y')
ylabel('time (ms), region X')


subplot(1,2,2)
imagesc(dcell.compsumGrayscale)
caxis([0 1])
axis square
titleuse=sprintf("post-inference \\sigma_{noise}=%d", sig);
title(titleuse, 'Interpreter', 'tex')
xlabel('time (ms), region Y')
ylabel('time (ms), region X')





