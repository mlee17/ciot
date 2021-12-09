% fitCumulativeGaussian
%
%      usage: fitCumulativeGaussian(x,y)
%         by: justin gardner
%       date: 10/11/2019
%    purpose: fit a cumulative gaussian to psychophysical data
%             y should be correct proportions (values between 0 and 1)
%    options: 'minParams=[-inf 0 0]' The minimum values parameters can go to
%             'maxParams=[inf inf inf] The maximum values parameters can go to
%             'initParams=[]' Starting values of parameters - defaults to median and std of x values
%             'maxIter=inf' How many iterations to search for best parameters
%             'fitType=basic' Fits just mean and standard deviation. Set to 'lapse' to also fit a lapse rate
%             'guessRate=0' Sets the guess rate, if you have a 2AFC you should set to 0.5
%
%       e.g.: 
%             x = [1 5 10 20 30 40];
%             y = [0.02 0.1 0.15 0.3 0.7 0.9];
%             fit = fitCumulativeGaussian(x,y,'fitType=lapse')
%             clf; plot(x,y,'ko'); hold on; plot(fit.fitX,fit.fitY,'r-');
%             xlabel('stimulus value');ylabel('Correct (proportion)');
%             title(sprintf('Mean %f Std %f lambda %f',fit.mean,fit.std,fit.lambda));
%              
%
%
function [bestfit_full, bestfit_amp, F_obt] = compareGaussianModels(x,y_high,y_med,y_low,varargin)

% check arguments
if nargin < 4
  help compareGaussianModels;
end

% parse input arguments
% getArgs(varargin,{'minParams',[-inf 0 -inf -inf],'maxParams',[inf inf inf inf],'initParams',[],'maxIter=inf'});

% make sure we have a column vector
x = x(:)';y_high = y_high(:)'; y_med = y_med(:)'; y_low = y_low(:)';

% set the initial parameters
% Full model: mean_high, mean_med, mean_low, std_high, std_med, std_low,
% amp_high, amp_med, amp_low, offset_high, offset_med, offset_low
minParams_full = [-inf -inf -inf 0 0 0 -inf -inf -inf -inf -inf -inf];
maxParams_full = [inf inf inf inf inf inf inf inf inf inf inf inf];
initParams_full = [median(x) median(x) median(x) 20 20 20 1 1 1 0 0 0]; %median(x) std(x) amplitude offset
nParams_full = 12; % 4*3
% Reduced (amplitude) model: mean_high, mean_med, mean_low, std, 
% amp_high, amp_med, amp_low, offset_high, offset_med, offset_low 
minParams_amp = [-inf -inf -inf 0 -inf -inf -inf -inf -inf -inf];
maxParams_amp = [inf inf inf inf inf inf inf inf inf inf];
initParams_amp = [median(x) median(x) median(x) 20 1 1 1 0 0 0];
nParams_amp = 10;

nObs = length(y_high) + length(y_med) + length(y_low);
  
% set optimization parametrs
optimParams = optimset('MaxIter',maxIter, 'TolFun',1e-8);

% <<< Full Model >>>
% some globals to keep track of what lsqnonlin does
global numIters;numIters = 0;
% fit function using lsqnonlin in LevenbergMarquardt mode.
[fitParams resnorm residual exitflag output lambda jacobian] = lsqnonlin(@gaussianErr_Full,initParams_full,minParams_full,maxParams_full,optimParams,x,y_high,y_med,y_low);
% [bestfit.fitparams, bestfit.fval, bestfit.exitflag] = fminsearch(@cgaussianErr, initParams, optimParams,x,y);

% Taken from Numerical Recipies, 
% the leastsq function seems to return the transposed gradient
% instead of the jacobian...
jacobian = jacobian'*jacobian;
reducedChiSquared = (residual*residual')/(length(y)-length(initParams));
covar = sqrt(reducedChiSquared * inv(jacobian));

% if we have the best fit then keep it.
bestfit_full = extractParams(fitParams);
bestfit_full.params = fitParams;
bestfit_full.covar = covar;
bestfit_full.output = output;

% get the fit
[bestfit_full.err bestfit_full.fit_high bestfit_full.fit_med bestfit_full.fit_low] = gaussianErr_Full(bestfit_full.params,x,y_high,y_med,y_low);

% compute r2 of fit
bestfit_full.r2 = 1-var(bestfit_full.err)/var(y);

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
jacobian = jacobian'*jacobian;
reducedChiSquared = (residual*residual')/(length(y)-length(initParams));
covar = sqrt(reducedChiSquared * inv(jacobian));

% if we have the best fit then keep it.
bestfit_amp = extractParams(fitParams);
bestfit_amp.params = fitParams;
bestfit_amp.covar = covar;
bestfit_amp.output = output;

% get the fit
[bestfit_amp.err bestfit_amp.fit_high bestfit_amp.fit_med bestfit_amp.fit_low] = gaussianErr_Amplitude(bestfit_amp.params,x,y_high,y_med,y_low);

% compute r2 of fit
bestfit_amp.r2 = 1-var(bestfit_amp.err)/var(y);

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
F_obt = ((SSE_amp - SSE_full) / (nParams_full - nParams_amp)) / (SSE_full / (nObs - nParams_full));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    GaussianErr    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [err, fit_high, fit_med, fit_low] = gaussianErr_Full(fitParams,x,y_high,y_med,y_low)

% get the parmas
p = extractParams(fitParams, 'Full');

% calculate the gaussian
% fit = (1/(p.std*sqrt(2*pi)))*exp(-1/2 * ((x-p.mean)/p.std).^2);
fit_high = p.amp_high * normpdf(x, p.mean_high, p.std_high) + p.offset_high;
fit_med = p.amp_med * normpdf(x, p.mean_med, p.std_med) + p.offset_med;
fit_low = p.amp_low * normpdf(x, p.mean_low, p.std_low) + p.offset_low;
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
fit_high = p.amp_high * normpdf(x, p.mean_high, p.std) + p.offset_high;
fit_med = p.amp_med * normpdf(x, p.mean_med, p.std) + p.offset_med;
fit_low = p.amp_low * normpdf(x, p.mean_low, p.std) + p.offset_low;
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
    p.mean_high = fitParams(1);
    p.mean_med = fitParams(2);
    p.mean_low = fitParams(3);
    p.std_high = fitParams(4);
    p.std_med = fitParams(5);
    p.std_low = fitParams(6);
    p.amp_high = fitParams(7);
    p.amp_med = fitParams(8);
    p.amp_low = fitParams(9);
    p.offset_high = fitParams(10);
    p.offset_med = fitParams(11);
    p.offset_low = fitParams(12);
elseif strcmp(fitType, 'Amplitude')
    p.mean_high = fitParams(1);
    p.mean_med = fitParams(2);
    p.mean_low = fitParams(3);
    p.std = fitParams(4);
    p.amp_high = fitParams(5);
    p.amp_med = fitParams(6);
    p.amp_low = fitParams(7);
    p.offset_high = fitParams(8);
    p.offset_med = fitParams(9);
    p.offset_low = fitParams(10);
else
  disp(sprintf('(compareGaussianModels:extractParams) Unknown fitType: %s',fitType));
  keyboard
end
      
      
