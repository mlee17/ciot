
%
function [bestfit_full, bestfit_amp, F_obt, Fr, AIC_full, AIC_amp] = compareGaussianModels_fixedOffset(x,y_high,y_med,y_low,varargin)

% check arguments
if nargin < 4
  help compareGaussianModels;
end

% parse input arguments
% getArgs(varargin,{'minParams',[-inf 0 -inf -inf],'maxParams',[inf inf inf inf],'initParams',[],'maxIter=inf'});

% make sure we have a column vector
x = x(:)';y_high = y_high(:)'; y_med = y_med(:)'; y_low = y_low(:)';

% set the initial parameters
% Full model:  std_high, std_med, std_low,
% amp_high, amp_med, amp_low, offset
minParams_full = [ 0 0 0 -inf -inf -inf -inf];
maxParams_full = [inf inf inf inf inf inf inf];
initParams_full = [20 20 20 1 1 1 0]; %median(x) std(x) amplitude offset
nParams_full = 7; % 4*3 -2
% Reduced (amplitude) model: mean_high, mean_med, mean_low, std, 
% amp_high, amp_med, amp_low, offset
minParams_amp = [0 -inf -inf -inf -inf];
maxParams_amp = [inf inf inf inf inf];
initParams_amp = [20 1 1 1 0];
nParams_amp = 5;

nObs = length(y_high) + length(y_med) + length(y_low);
  
% set optimization parametrs
maxIter = inf;
optimParams = optimset('MaxIter',maxIter, 'TolFun',1e-8, 'MaxFunEvals', inf);

% <<< Full Model >>>
% some globals to keep track of what lsqnonlin does
global numIters;numIters = 0;
% fit function using lsqnonlin in LevenbergMarquardt mode.
[fitParams resnorm residual exitflag output lambda jacobian] = lsqnonlin(@gaussianErr_Full,initParams_full,minParams_full,maxParams_full,optimParams,x,y_high,y_med,y_low);
% [bestfit.fitparams, bestfit.fval, bestfit.exitflag] = fminsearch(@cgaussianErr, initParams, optimParams,x,y);

% Taken from Numerical Recipies, 
% the leastsq function seems to return the transposed gradient
% instead of the jacobian...
% jacobian = jacobian'*jacobian;
% reducedChiSquared = (residual*residual')/(length(y)-length(initParams));
% covar = sqrt(reducedChiSquared * inv(jacobian));

% if we have the best fit then keep it.
bestfit_full = extractParams(fitParams,'Full');
bestfit_full.params = fitParams;
% bestfit_full.covar = covar;
bestfit_full.output = output;

% get the fit
[bestfit_full.err bestfit_full.fit_high bestfit_full.fit_med bestfit_full.fit_low] = gaussianErr_Full(bestfit_full.params,x,y_high,y_med,y_low);

% compute r2 of fit
bestfit_full.r2 = 1-var(bestfit_full.err)/(var([y_high; y_med; y_low]));

SSE_full = 0;
for i = 1:length(bestfit_full.err)
    SSE_full = SSE_full + bestfit_full.err(i)^2;
end

% compute a smoother fit (i.e. nFitPoints along x axis)
nFitPoints = 1000;
bestfit_full.fitX = min(x):(max(x)-min(x))/(nFitPoints-1):max(x);

% note here that the y variable is just a dummy value since we don't
% care about the error 
[~,bestfit_full.fitY_high,bestfit_full.fitY_med,bestfit_full.fitY_low] = gaussianErr_Full(bestfit_full.params,bestfit_full.fitX,zeros(1,nFitPoints),zeros(1,nFitPoints),zeros(1,nFitPoints));

% <<< Amplitude (Reduced) Model >>>
% some globals to keep track of what lsqnonlin does
numIters = 0;
% fit function using lsqnonlin in LevenbergMarquardt mode.
[fitParams resnorm residual exitflag output lambda jacobian] = lsqnonlin(@gaussianErr_Amplitude,initParams_amp,minParams_amp,maxParams_amp,optimParams,x,y_high,y_med,y_low);
% [bestfit.fitparams, bestfit.fval, bestfit.exitflag] = fminsearch(@cgaussianErr, initParams, optimParams,x,y);

% Taken from Numerical Recipies, 
% the leastsq function seems to return the transposed gradient
% instead of the jacobian...
% jacobian = jacobian'*jacobian;
% reducedChiSquared = (residual*residual')/(length(y)-length(initParams));
% covar = sqrt(reducedChiSquared * inv(jacobian));

% if we have the best fit then keep it.
bestfit_amp = extractParams(fitParams, 'Amplitude');
bestfit_amp.params = fitParams;
% bestfit_amp.covar = covar;
bestfit_amp.output = output;

% get the fit
[bestfit_amp.err bestfit_amp.fit_high bestfit_amp.fit_med bestfit_amp.fit_low] = gaussianErr_Amplitude(bestfit_amp.params,x,y_high,y_med,y_low);

% compute r2 of fit
bestfit_amp.r2 = 1-var(bestfit_amp.err)/(var([y_high; y_med; y_low]));

SSE_amp = 0;
for i = 1:length(bestfit_amp.err)
    SSE_amp = SSE_amp + bestfit_amp.err(i)^2;
end

% compute a smoother fit (i.e. nFitPoints along x axis)
nFitPoints = 1000;
bestfit_amp.fitX = min(x):(max(x)-min(x))/(nFitPoints-1):max(x);

% note here that the y variable is just a dummy value since we don't
% care about the error 
[~,bestfit_amp.fitY_high,bestfit_amp.fitY_med,bestfit_amp.fitY_low] = gaussianErr_Amplitude(bestfit_amp.params,bestfit_amp.fitX,zeros(1,nFitPoints),zeros(1,nFitPoints),zeros(1,nFitPoints));

%%%%%%%%%
% Compute F Statistics
F_obt = ((SSE_amp - SSE_full) / (nParams_full - nParams_amp)) / (SSE_full / (nObs - nParams_full - 1));
Fr = ((bestfit_full.r2 - bestfit_amp.r2) / (nParams_full - nParams_amp)) / ((1-bestfit_full.r2) / (nObs - nParams_full - 1));
bestfit_full.Fr = Fr;
bestfit_amp.Fr = Fr;

bestfit_amp.F_obt = F_obt;
df1 = nParams_full - nParams_amp;
df2 = nObs - nParams_full;
bestfit_amp.df1 = df1;
bestfit_amp.df2 = df2;

% compute AIC
% AIC ~= 2k + n*ln(RSS) % k = n params, n = n observations
AIC_full = 2*nParams_full + nObs*log(SSE_full);
AIC_amp = 2*nParams_amp + nObs*log(SSE_amp);
bestfit_full.AIC = AIC_full;
bestfit_amp.AIC = AIC_amp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    GaussianErr    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [err, fit_high, fit_med, fit_low] = gaussianErr_Full(fitParams,x,y_high,y_med,y_low)

% get the parmas
p = extractParams(fitParams, 'Full');

% calculate the gaussian
% fit = (1/(p.std*sqrt(2*pi)))*exp(-1/2 * ((x-p.mean)/p.std).^2);
fit_high = p.amp_high * normpdf(x, 0, p.std_high) + p.offset;
fit_med = p.amp_med * normpdf(x, 0, p.std_med) + p.offset;
fit_low = p.amp_low * normpdf(x, 0, p.std_low) + p.offset;
% normpdf(x,p.mean,p.std)

% update number of iterations
global numIters;
numIters = numIters+1;

err = [y_high-fit_high; y_med-fit_med; y_low-fit_low];

function [err, fit_high, fit_med, fit_low] = gaussianErr_Amplitude(fitParams,x,y_high,y_med,y_low)

% get the parmas
p = extractParams(fitParams, 'Amplitude');

% calculate the gaussian
% fit = (1/(p.std*sqrt(2*pi)))*exp(-1/2 * ((x-p.mean)/p.std).^2);
fit_high = p.amp_high * normpdf(x, 0, p.std) + p.offset;
fit_med = p.amp_med * normpdf(x, 0, p.std) + p.offset;
fit_low = p.amp_low * normpdf(x, 0, p.std) + p.offset;
% normpdf(x,p.mean,p.std)

% update number of iterations
global numIters;
numIters = numIters+1;

err = [y_high-fit_high; y_med-fit_med; y_low-fit_low];

%%%%%%%%%%%%%%%%%%%%%%%
%    extractParams    %
%%%%%%%%%%%%%%%%%%%%%%%
function p = extractParams(fitParams, fitType)

% extrat the parameters
  
if strcmp(fitType, 'Full')
   
    p.std_high = fitParams(1);
    p.std_med = fitParams(2);
    p.std_low = fitParams(3);
    p.amp_high = fitParams(4);
    p.amp_med = fitParams(5);
    p.amp_low = fitParams(6);
    p.offset = fitParams(7);
    
elseif strcmp(fitType, 'Amplitude')
   
    p.std = fitParams(1);
    p.amp_high = fitParams(2);
    p.amp_med = fitParams(3);
    p.amp_low = fitParams(4);
    p.offset = fitParams(5);
    
else
  disp(sprintf('(compareGaussianModels:extractParams) Unknown fitType: %s',fitType));
  keyboard
end
      
      
