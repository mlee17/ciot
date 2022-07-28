function simulateTuning
%% make a gabor receptive field
rf.size = 2;
rf.sd = 0.5;
rf.sf = 1;

orientations = [0 7 15 30 45 75 90] + 90;
rf.tuning = [0 30 60 90 120 150] + 90;

contrasts = [0.3 0.15 0.075];
mglOpen;
grating = mglMakeGrating(rf.size,rf.size, rf.sf, rf.tuning(1), 0);
gaussian = mglMakeGaussian(rf.size,rf.size, rf.sd, rf.sd);
defaultRF = grating.*gaussian;
for i = 1:length(rf.tuning)
    tempRF = mglMakeGrating(rf.size, rf.size, rf.sf, rf.tuning(i), 0);
    receptiveField{i} = tempRF.*gaussian;
end
for i = 1:length(orientations)
    tempStim = mglMakeGrating(rf.size,rf.size, rf.sf, orientations(i), 0);
    stim{i} = tempStim.*gaussian;
end
mglClose;
mult = defaultRF.*stim{1}*contrasts(1);
maxResponse = sum(mult(:));
for c = 1:length(contrasts)
for i = 1:length(orientations)
    mult = defaultRF.*stim{i}*contrasts(c);
    amplitudeMask{c}(i) = sum(mult(:)) / maxResponse;
end
    amplitudeTest(c) = amplitudeMask{c}(1);
end


%% input signals
deltaT = 0.001;
t = 0:deltaT:1;%time
f1 = 3;
f2 = 5;
f1wave = sin(2*pi*f1*t);
f2wave = sin(2*pi*f2*t);

figure; hold on;
plot(t, f1wave);
plot(t, f2wave);

f1waveFullRect = f1wave;%abs(f1wave); % full rect
% f1wave(f1wave<0) = 0; %half rect
f2waveFullRect = f2wave;%abs(f2wave);

% abs(f1wave+f2wave);

%% normalization
sigma = 0.2; % semi-saturation paramter
p = 2;% power for numerator
q = 2; % power for denominator

for c = 1:length(contrasts)
for s = 1:length(orientations) % for each orientation offset
    
numerator{c}{s} = (amplitudeTest(c)*f1wave + amplitudeMask{c}(s)*f2wave).^p;%abs(amplitudeTest(c)*f1wave + amplitudeMask{c}(s)*f2wave).^p;
normalizationPool{c}{s} = zeros(1,length(t));
for i = 1:length(rf.tuning) % for each receptiveField (normalization pool)
    mult_test = receptiveField{i} .* stim{1} *contrasts(c); % fixed to stim 1 (test)
    thisAmp_test = sum(mult_test(:)) / maxResponse;
    mult_mask = receptiveField{i} .* stim{s} *contrasts(c);
    thisAmp_mask = sum(mult_mask(:)) / maxResponse;
    
%     if orientations(s) > 30 + 90
        thisNorm = (thisAmp_test*f1wave).^q + (thisAmp_mask*f2wave).^q;%abs(thisAmp_test*f1wave).^q + abs(thisAmp_mask*f2wave).^q;
%     else
%         thisNorm = abs(thisAmp_test*f1wave + thisAmp_mask*f2wave).^q;
%     end
    
    normalizationPool{c}{s} = normalizationPool{c}{s} + thisNorm;
end
Response{c}{s} = numerator{c}{s} ./ (normalizationPool{c}{s} + sigma^q);
% figure;
fftMag{c}{s} = 2*deltaT*abs(fft(Response{c}{s}));
% stem(1:50, fftMag(2:51));

secondDiff{c}(s) = fftMag{c}{s}(3);
secondSum{c}(s) = fftMag{c}{s}(9);

fftMag_numerator{c}{s} = 2*deltaT*abs(fft(numerator{c}{s}));
secondDiff2{c}(s) = fftMag_numerator{c}{s}(3);
secondSum2{c}(s) = fftMag_numerator{c}{s}(9);

end
end
%%
figure;
stem(1:40, fftMag{1}{1}(2:41), 'Color', [0.5 0.5 0.5],'LineWidth',1.5);
hold on;
stem([6,10], [fftMag{1}{1}(7) fftMag{1}{1}(11)],'filled', 'Color', [0.2 0.3 0.8],'LineWidth',3);
stem([2,8], [fftMag{1}{1}(3) fftMag{1}{1}(9)], 'filled','Color', [0.75 0.3 0.55],'LineWidth',3);
stem([4,16], [fftMag{1}{1}(5) fftMag{1}{1}(17)], 'filled','Color',[0.9 0.7 0.8],'LineWidth',3);
box off;
xlabel('Frequency (Hz)');
ylabel('Amplitude');
set(gca,'ytick', []);

%%
for c = 1:length(contrasts)
temp = flip(secondDiff{c});
secondDiffAll{c} = [temp, secondDiff{c}(2:end)];
temp = flip(secondSum{c});
secondSumAll{c} = [temp, secondSum{c}(2:end)];
end
temp = flip(-(orientations-90));
orientationAll = [temp, orientations(2:end)-90];

colors = brewermap(7, 'BuPu');

figure;hold on;
plot(orientationAll, secondSumAll{1},'-o', 'Color',colors(7,:), 'LineWidth',3); %[0.75 0.3 0.55]
plot(orientationAll, secondSumAll{2},'-o', 'Color',colors(5,:), 'LineWidth',3); %[0.75 0.3 0.55]
plot(orientationAll, secondSumAll{3},'-o', 'Color',colors(3,:), 'LineWidth',3); %[0.75 0.3 0.55]
myerrorbar(orientationAll, secondSumAll{1}, 'Color', colors(7,:), 'Symbol','o', 'MarkerSize', 10);
myerrorbar(orientationAll, secondSumAll{2}, 'Color', colors(5,:), 'Symbol','o', 'MarkerSize', 10);
myerrorbar(orientationAll, secondSumAll{3}, 'Color', colors(3,:), 'Symbol','o', 'MarkerSize', 10);

xlabel('Orienation Offset (Deg)');
ylabel('f1+f2 Amplitude');
drawPublishAxis;

figure;hold on;
plot(orientationAll, secondDiffAll{1},'-o', 'Color',colors(7,:), 'LineWidth',3); %[0.75 0.3 0.55]
plot(orientationAll, secondDiffAll{2},'-o', 'Color',colors(5,:), 'LineWidth',3); %[0.75 0.3 0.55]
plot(orientationAll, secondDiffAll{3},'-o', 'Color',colors(3,:), 'LineWidth',3); %[0.75 0.3 0.55]
myerrorbar(orientationAll, secondDiffAll{1}, 'Color', colors(7,:), 'Symbol','o', 'MarkerSize', 10);
myerrorbar(orientationAll, secondDiffAll{2}, 'Color', colors(5,:), 'Symbol','o', 'MarkerSize', 10);
myerrorbar(orientationAll, secondDiffAll{3}, 'Color', colors(3,:), 'Symbol','o', 'MarkerSize', 10);
xlabel('Orienation Offset (Deg)');
ylabel('f2-f1 Amplitude');
drawPublishAxis;


%%
% f1 + f2 low contrast std = 0.09
% f2 - f1 low contrast std = 0.17
std1 = 0.09;
std2 = 0.17;
nSubj = [5, 10, 20, 30];
for j = 1:length(nSubj)
    for iSubj = 1:nSubj(j)
        for i = 1:length(secondDiffAll{3})
            subjData{j}(iSubj,i) = secondDiffAll{3}(i) + std2*randn(1);  
        end
    end
    meanData{j} = mean(subjData{j},1);
    stdData{j} = std(subjData{j},1);
end

nColors = brewermap(8, 'OrRd');

figure; hold on;
plot(orientationAll, secondDiffAll{3},'-', 'Color',colors(3,:), 'LineWidth',3); %[0.75 0.3 0.55]
myerrorbar(orientationAll, secondDiffAll{3}, 'Color', colors(3,:), 'Symbol','o', 'MarkerSize', 10);

for j = 1:length(nSubj)
    myerrorbar(orientationAll, meanData{j}, 'yError', stdData{j}, 'Color', nColors(j+4,:), 'Symbol','--','MarkerSize',10);
end





