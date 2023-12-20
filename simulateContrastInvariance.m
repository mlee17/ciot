function simulateContrastInvariance
% parameters for the simple cell
meshsize = 64;
sfPref = 4;
orientationPref = 90;
c50 = .1;

% make a simple with normalization
s = simpleCellWithNormalization(meshsize, orientationPref, sfPref, c50);

% init response array
r = [];

% make a grating that has the preferred orientation and spatial frequency
% of the model cell
stimulus = makeGrating(s.x, s.y, orientationPref, 0, sfPref);

% loop over contrasts
nContrasts = 500;
contrasts = linspace(0,1,nContrasts);
for iContrast = 1:nContrasts
    r(iContrast) = computeSimpleCellWithNormalizationResponse(s, contrasts(iContrast)*stimulus);
end

% plot
figure;
semilogx(contrasts,r);
xlabel('Contrast');
ylabel('Response');


% contrast-invariant orientation tuning
% init response array
nContrasts = 3;
nOrientations = 180;
r = zeros(nContrasts,nOrientations);

% loop over contrasts
contrasts = [0.02 0.04 0.06];
for iContrast = 1:nContrasts
    % loop over orientations
    orientations = 1:1:180;
    for orientation = orientations
        % make a grating of the appropirate orientation
        stimulus = makeGrating(s.x, s.y, orientation, 0, sfPref);
        
        % compute simple cell with normalization response
        r(iContrast,orientation) = computeSimpleCellWithNormalizationResponse(s, contrasts(iContrast) * stimulus);
        
    end
end

colors = brewermap(7, 'BuPu');
figure;
plot(r(1,:),'Color',colors(3,:), 'LineWidth',3);
hold on;
plot(r(2,:),'Color',colors(5,:), 'LineWidth',3);
plot(r(3,:),'Color',colors(7,:), 'LineWidth',3);
xlabel('Orientation (deg)');
ylabel('Neural Response (spikes/sec)');

 % contrast orientation tuning without normalization
 r2 = zeros(nContrasts,nOrientations);

 % loop over contrasts
 for iContrast = 1:nContrasts
     % loop over orientations
     orientations = 1:1:180;
     for orientation = orientations
         % make a grating of the appropirate orientation
         stimulus = makeGrating(s.x, s.y, orientation, 0, sfPref);
         % compute simple cell without normalization response
         r2(iContrast, orientation) = computeSimpleCellWithoutNormalizationResponse(s, contrasts(iContrast)*stimulus);
     end
 end
 
 figure;
 plot(r2(1,:),'Color',colors(3,:), 'LineWidth',3);
hold on;
plot(r2(2,:),'Color',colors(5,:), 'LineWidth',3);
plot(r2(3,:),'Color',colors(7,:), 'LineWidth',3);
xlabel('Orientation (deg)');
ylabel('Neural Response (spikes/sec)');


keyboard

function [x,y,extents] = getMeshPoints(nPoints)
x = linspace(-1,1,nPoints);
y = linspace(-1,1,nPoints);
extents = [min(x), max(x), min(y), max(y)];
[x,y] = meshgrid(x,y);

function gaussian = makeGaussian(x,y,sigma)
gaussian = exp(-(x.^2+y.^2)/(2*sigma^2));

function grating = makeGrating(x,y, orientation, spatialPhase, spatialFrequency)
% we wil convert orientation and spatialPhase into radians
orientation = pi*orientation/180;
spatialPhase = pi*spatialPhase/180;
% we need to convert spatial frequency into cycles/image
% remember that we made the extents in getMeshPoints
% to go from -1 to 1, so we want that to go from -pi to pi
spatialFrequency = spatialFrequency * pi;

% make the grating
grating = cos(spatialFrequency*(x.*cos(orientation)+y.*sin(orientation))+spatialPhase);

function simple = simpleCell(meshsize, orientationPreference, spatialPhase, spatialFrequencyPreference, exponent)
simple.meshsize = meshsize;

[simple.x, simple.y, simple.extents] = getMeshPoints(meshsize);

% compute the grating needed
grating = makeGrating(simple.x, simple.y, orientationPreference, spatialPhase, spatialFrequencyPreference);

% compute the gaussian (we fix the size here, but of course that could be a
% passed in parameter
gaussian = makeGaussian(simple.x, simple.y, 0.2);

% now we can make the garbor receptive field
simple.linearRF = grating .* gaussian;

%the exponent is just something we store
simple.exponent = exponent;

function response = computeSimpleCellResponse(simple, stimulus)
%  well, the first step is to apply the linear receptive field
%  which means to take the dot product of the linear receptive field
%  and the stimulus. Note, I don't understand what numpy's dot product
%  does in 2D - so maybe there is a better way to write this. Instead
%  I'm going to go really basic here. Take the element-wise multiplication
%  of stimulus and RF and then add that all up together
response = simple.linearRF .* stimulus;
response = sum(response);

% since the units of the output are aribtrary, let's
% make them a bit more intepretable where 1 would be the
% maximum possible reponse of the RF with it's most preferred
% stimulus. What is the most preferred stimulus, well the one
% that exactly matches the RF!
maxResponse = simple.linearRF .* simple.linearRF;
maxResponse = sum(maxResponse);

% now normalize by this max response
response = response / maxResponse;

% now we apply a threshold
if response < 0
    response = 0;
end

% and apply the static non-linearity
response = response^simple.exponent;

function simpleNorm = simpleCellWithNormalization(meshsize,orientationPreference, spatialFrequencyPreference, c50)
simpleNorm.meshsize = meshsize;
[simpleNorm.x, simpleNorm.y, simpleNorm.extents] = getMeshPoints(meshsize);

% create the simple cell receptive field - note that we are leaving out the exponent
% because we are going to compute that ourselves (i.e. setting exponent to 1)
simpleNorm.rf = simpleCell(meshsize,orientationPreference,0,spatialFrequencyPreference,1);

% now create simple cells for the normalization pool
simpleNorm.nNormPool = 12;
% simpleNorm.normPool = [];
for i = 1:simpleNorm.nNormPool
    % get an orientation preference
    normOrientationPreference = i*360/simpleNorm.nNormPool;
    % create the appropriate simple cell
    simpleNorm.normPool(i) = simpleCell(meshsize,normOrientationPreference,0,spatialFrequencyPreference,1);
end

% here, just trying to set the range of the norm pool similarly to the simple cell response
% maximum resonse of 1 for the best possible stimulus, so to do that, we use the same trick we have done
% with a simple cell. We compute the best possible stimulus for one of the
% filters (it's matched stimulus) and then hit all the other receptive fields
% with the same stimulus and see what the summed response is. This summed
% value should be the highest possible response you can get with the norm
% pool. nb I'm not quite sure I have all of this right. It is possible
% for the stimulus to give a larger response than the maxResponse - i.e. if you take the
% preferred stimulus and up the contrast. I think the most important thing is just
% that the response from the normalization pool is similar to the linear receptive field
% as a function of the contrast, so that this behaves like the naka-rushton eqution
% (i.e. if the normpool grows 10 times as fast as a function of contrast then the
% linear receptive fields, it wouldn't be implementing c^2 / c50^2 +c^2), but would
% be implementing c^2 / c50^2 to (10c)^2 . M
% maxResponse = [];
for i = 1:simpleNorm.nNormPool
    maxResponse(i) = computeSimpleCellResponse(simpleNorm.normPool(i), simpleNorm.normPool(1).linearRF);
end
simpleNorm.normPoolMaxResponse = sum(maxResponse);

% keep the c50 value
simpleNorm.c50 = c50;

% set the exponent to be 2
simpleNorm.exponent = 2;


function response = computeSimpleCellWithNormalizationResponse(s, stimulus)
% compute response to central receptive field
rfResponse = computeSimpleCellResponse(s.rf, stimulus);

% compute normalization pool response
normPoolResponse = [];
for i = 1:s.nNormPool
    normPoolResponse(i) = computeSimpleCellResponse(s.normPool(i), stimulus);
end

% sum the norm pool
normPoolResponse = sum(normPoolResponse);

% and divide by the maximum response - just to get the output of the
% normalization pool into the same range as the linear receptive field
% see note above
normPoolResponse = normPoolResponse / s.normPoolMaxResponse;

% compute normalized response - there it is - the Naka-Rushton equation that
% makes this whole thing work!
response = rfResponse.^s.exponent./(s.c50.^s.exponent + normPoolResponse.^s.exponent);

function response = computeSimpleCellWithoutNormalizationResponse(s, stimulus)
% compute response to central receptive field
rfResponse = computeSimpleCellResponse(s.rf, stimulus);

% compute response which just maxes out at some maximum
responseMax = 0.05;
if rfResponse > responseMax
    response = responseMax;
else
    response = rfResponse;
end





