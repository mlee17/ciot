
%
function [bestfit_full, bestfit_amp, F_obt, AIC_full, AIC_amp] = compareGaussianModels_IM(x,y_high,y_med,y_low, y_high2, y_med2, y_low2, varargin)

% check arguments
if nargin < 4
  help compareGaussianModels;
end

% parse input arguments
% getArgs(varargin,{'minParams',[-inf 0 -inf -inf],'maxParams',[inf inf inf inf],'initParams',[],'maxIter=inf'});

% make sure we have a column vector
x = x(:)';y_high = y_high(:)'; y_med = y_med(:)'; y_low = y_low(:)';
y_high2 = y_high2(:)'; y_med2 = y_med2(:)'; y_low2 = y_low2(:)';
% set the initial parameters
% Full model: std_high1, std_med1, std_low1, std_high2, std_med2, std_low2
% amp_high1, amp_med1, amp_low1, amp_high2, amp_med2, amp_low2 offset1, offset2
minParams_full = [0 0 0 0 0 0 -inf -inf -inf -inf -inf -inf -inf -inf];
maxParams_full = [inf inf inf inf inf inf inf inf inf inf inf inf inf inf];
initParams_full = [20 20 20 20 20 20 1 1 1 1 1 1 0 0]; % std(x) amplitude offset
nParams_full = length(minParams_full); 
% Reduced (amplitude) model: std, 
% amp_high, amp_med, amp_low, amp_high2, amp_med2, amp_low2, 
% offset1, offset2
minParams_amp = [0 -inf -inf -inf -inf -inf -inf -inf -inf];
maxParams_amp = [inf inf inf inf inf inf inf inf inf];
initParams_amp = [20 1 1 1 1 1 1 0 0];
nParams_amp = length(minParams_amp);

nObs = length(y_high) + length(y_med) + length(y_low) + length(y_high2) + length(y_med2) + length(y_low2);
  
% set optimization parametrs
maxIter = inf;
optimParams = optimset('MaxIter',maxIter, 'TolFun',1e-8,'MaxFunEvals',inf);

% <<< Full Model >>>
% some globals to keep track of what lsqnonlin does
global numIters;numIters = 0;
% fit function using lsqnonlin in LevenbergMarquardt mode.
[fitParams resnorm residual exitflag output lambda jacobian] = lsqnonlin(@gaussianErr_Full,initParams_full,minParams_full,maxParams_full,optimParams,x,y_high,y_med,y_low, y_high2,y_med2,y_low2);
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
[bestfit_full.err bestfit_full.fit_high bestfit_full.fit_med bestfit_full.fit_low bestfit_full.fit_high2 bestfit_full.fit_med2 bestfit_full.fit_low2] = ...
    gaussianErr_Full(bestfit_full.params,x,y_high,y_med,y_low, y_high2,y_med2,y_low2);

% compute r2 of fit
% bestfit_full.r2 = 1-var(bestfit_full.err)/var(y);

SSE_full = 0;
for i = 1:length(bestfit_full.err)
    SSE_full = SSE_full + bestfit_full.err(i)^2;
end

% compute a smoother fit (i.e. nFitPoints along x axis)
nFitPoints = 1000;
bestfit_full.fitX = min(x):(max(x)-min(x))/(nFitPoints-1):max(x);

% note here that the y variable is just a dummy value since we don't
% care about the error 
[~,bestfit_full.fitY_high,bestfit_full.fitY_med,bestfit_full.fitY_low, bestfit_full.fitY_high2,bestfit_full.fitY_med2,bestfit_full.fitY_low2] =...
    gaussianErr_Full(bestfit_full.params,bestfit_full.fitX,zeros(1,nFitPoints),zeros(1,nFitPoints),zeros(1,nFitPoints),zeros(1,nFitPoints),zeros(1,nFitPoints),zeros(1,nFitPoints));

% <<< Amplitude (Reduced) Model >>>
% some globals to keep track of what lsqnonlin does
numIters = 0;
% fit function using lsqnonlin in LevenbergMarquardt mode.
[fitParams resnorm residual exitflag output lambda jacobian] = lsqnonlin(@gaussianErr_Amplitude,initParams_amp,minParams_amp,maxParams_amp,optimParams,x,y_high,y_med,y_low, y_high2,y_med2,y_low2);
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
[bestfit_amp.err bestfit_amp.fit_high bestfit_amp.fit_med bestfit_amp.fit_low bestfit_amp.fit_high2 bestfit_amp.fit_med2 bestfit_amp.fit_low2] = ...
    gaussianErr_Amplitude(bestfit_amp.params,x,y_high,y_med,y_low, y_high2,y_med2,y_low2);

% compute r2 of fit
% bestfit_amp.r2 = 1-var(bestfit_amp.err)/var(y);

SSE_amp = 0;
for i = 1:length(bestfit_amp.err)
    SSE_amp = SSE_amp + bestfit_amp.err(i)^2;
end

% compute a smoother fit (i.e. nFitPoints along x axis)
nFitPoints = 1000;
bestfit_amp.fitX = min(x):(max(x)-min(x))/(nFitPoints-1):max(x);

% note here that the y variable is just a dummy value since we don't
% care about the error 
[~,bestfit_amp.fitY_high,bestfit_amp.fitY_med,bestfit_amp.fitY_low, bestfit_amp.fitY_high2,bestfit_amp.fitY_med2,bestfit_amp.fitY_low2] =...
    gaussianErr_Amplitude(bestfit_amp.params,bestfit_amp.fitX,zeros(1,nFitPoints),zeros(1,nFitPoints),zeros(1,nFitPoints),zeros(1,nFitPoints),zeros(1,nFitPoints),zeros(1,nFitPoints));

%%%%%%%%%
% Compute F Statistics
F_obt = ((SSE_amp - SSE_full) / (nParams_full - nParams_amp)) / (SSE_full / (nObs - nParams_full - 1));

% compute AIC
AIC_full = 2*nParams_full + nObs*log(SSE_full);
AIC_amp = 2*nParams_amp + nObs*log(SSE_amp);
bestfit_full.AIC = AIC_full;
bestfit_amp.AIC = AIC_amp;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    GaussianErr    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [err, fit_high, fit_med, fit_low, fit_high2, fit_med2, fit_low2] = gaussianErr_Full(fitParams,x,y_high,y_med,y_low,y_high2,y_med2,y_low2)

% get the parmas
p = extractParams(fitParams, 'Full');

% calculate the gaussian
% fit = (1/(p.std*sqrt(2*pi)))*exp(-1/2 * ((x-p.mean)/p.std).^2);
fit_high = p.amp_high * normpdf(x, 0, p.std_high) + p.offset;
fit_med = p.amp_med * normpdf(x, 0, p.std_med) + p.offset;
fit_low = p.amp_low * normpdf(x, 0, p.std_low) + p.offset;
fit_high2 = p.amp_high2 * normpdf(x, 0, p.std_high2) + p.offset2;
fit_med2 = p.amp_med2 * normpdf(x, 0, p.std_med2) + p.offset2;
fit_low2 = p.amp_low2 * normpdf(x, 0, p.std_low2) + p.offset2;
% normpdf(x,p.mean,p.std)

% update number of iterations
global numIters;
numIters = numIters+1;

err = [y_high-fit_high; y_med-fit_med; y_low-fit_low; y_high2-fit_high2; y_med2-fit_med2; y_low2-fit_low2];

function [err, fit_high, fit_med, fit_low, fit_high2, fit_med2, fit_low2] = gaussianErr_Amplitude(fitParams,x,y_high,y_med,y_low,y_high2,y_med2,y_low2)

% get the parmas
p = extractParams(fitParams, 'Amplitude');

% calculate the gaussian
% fit = (1/(p.std*sqrt(2*pi)))*exp(-1/2 * ((x-p.mean)/p.std).^2);
fit_high = p.amp_high * normpdf(x, 0, p.std) + p.offset;
fit_med = p.amp_med * normpdf(x, 0, p.std) + p.offset;
fit_low = p.amp_low * normpdf(x, 0, p.std) + p.offset;
fit_high2 = p.amp_high2 * normpdf(x, 0, p.std) + p.offset2;
fit_med2 = p.amp_med2 * normpdf(x, 0, p.std) + p.offset2;
fit_low2 = p.amp_low2 * normpdf(x, 0, p.std) + p.offset2;
% normpdf(x,p.mean,p.std)

% update number of iterations
global numIters;
numIters = numIters+1;

err = [y_high-fit_high; y_med-fit_med; y_low-fit_low; y_high2-fit_high2; y_med2-fit_med2; y_low2-fit_low2];

%%%%%%%%%%%%%%%%%%%%%%%
%    extractParams    %
%%%%%%%%%%%%%%%%%%%%%%%
function p = extractParams(fitParams, fitType)

% extrat the parameters
  
if strcmp(fitType, 'Full')
    
    p.std_high = fitParams(1);
    p.std_med = fitParams(2);
    p.std_low = fitParams(3);
    p.std_high2 = fitParams(4);
    p.std_med2 = fitParams(5);
    p.std_low2 = fitParams(6);
    p.amp_high = fitParams(7);
    p.amp_med = fitParams(8);
    p.amp_low = fitParams(9);
    p.amp_high2 = fitParams(10);
    p.amp_med2 = fitParams(11);
    p.amp_low2 = fitParams(12);
    p.offset = fitParams(13);
    p.offset2 = fitParams(14);
    
elseif strcmp(fitType, 'Amplitude')
    
    p.std = fitParams(1);
    p.amp_high = fitParams(2);
    p.amp_med = fitParams(3);
    p.amp_low = fitParams(4);
    p.amp_high2 = fitParams(5);
    p.amp_med2 = fitParams(6);
    p.amp_low2 = fitParams(7);
    p.offset = fitParams(8);
    p.offset2 = fitParams(9);
    
else
  disp(sprintf('(compareGaussianModels:extractParams) Unknown fitType: %s',fitType));
  keyboard
end
      
      
