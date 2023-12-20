% coherent averaging

subjectList = {'nl-1624','nl-jg','nl-1928','nl-0055','nl-1909','nl-1908','nl-0051','nl-0034','nl-myl','nl-rta','nl-0057','nl-2274','nl-2275',...
    'nl-2277','nl-2278','nl-2279','nl-2280','nl-2281','nl-2282', 'nl-2283','nl-2285', 'nl-2287', 'nl-2300','nl-1539','nl-2113','nl-2126','nl-2215',...
    'nl-2329','nl-2331','nl-2332','nl-2286'};%,'nl-2335','nl-2336','nl-2343','nl-2344','nl-2346','nl-2347','nl-2349','nl-2350','nl-2351','nl-2353', 'nl-2374',...
    %'nl-2378', 'nl-2379','nl-2381', 'nl-2383'};
fileDir = '~/data/orientationRCA/snr';
dataDir = '~/data/orientationtuning';
matfileDir = 'Exp_MATL_HCN_128_Avg';
taskList = [zeros(31,1)];%;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1];

[orderedSNR, rankInd] = sort(meanSNR(1:nSubj,5), 'descend');
isOver2 = orderedSNR>2;

subjectList = subjectList(rankInd);
subjectList = subjectList(isOver2);
taskList = task(rankInd);
taskList = taskList(isOver2);

nSubj= length(subjectList);

target_frequencies = [8, 2];

output_cos = []; output_sin = []; output_amp=[]; output_side_cos=[]; output_side_sin=[];
for s = 1:length(subjectList)
    thisFolder = dir(fullfile(dataDir, [subjectList{s},'*']));
    thisDir = fullfile(dataDir, thisFolder.name, matfileDir);
        
    [output_cos{s}, output_sin{s}, output_amp{s}, output_side_cos{s}, output_side_sin{s}] =...
        getPhases_occipital(thisDir, target_frequencies, taskList(s));
    
end

%% Coherent Averaging

tempCos = []; %nFreq, nCond
tempSin = [];
sideCos = [];
sideSin = [];

for iFreq = 1:2
tempCosAll = [];
tempSinAll = [];
sideCosAll = [];
sideSinAll = [];
    for s = 1:length(subjectList)
        tempCosAll = [tempCosAll ; squeeze(output_cos{s}(iFreq,:))];
        tempSinAll = [tempSinAll ; squeeze(output_sin{s}(iFreq,:))];
        sideCosAll = [sideCosAll ; squeeze(output_side_cos{s}(iFreq,:))];
        sideSinAll = [sideSinAll ; squeeze(output_side_sin{s}(iFreq,:))];
    end
    tempCos(iFreq,:) = mean(tempCosAll,1);
    tempSin(iFreq,:) = mean(tempSinAll,1);
    sideCos(iFreq,:) = mean(sideCosAll,1);
    sideSin(iFreq,:) = mean(sideSinAll,1);
end

tempAmp = sqrt(tempCos.^2 + tempSin.^2);
tempAmpSide = sqrt(sideCos.^2 + sideSin.^2);

 freqLabels = {'1F1+1F2', '-1F1+1F2'};


HighContrast = [1:7];
MedContrast = [10:16];
LowContrast = [19:25];
ConditionsToAnalyze = [HighContrast,MedContrast,LowContrast];

orientationAll = [-90 -75 -45 -30 -15 -7 0 7 15 30 45 75 90];

highAmp = []; medAmp = []; lowAmp = [];
for iFreq = 1:2

p = 1;
  for i = HighContrast
    highAmp{iFreq}(p) = tempAmp(iFreq,i); % Row IM, Col RC
    highSide{iFreq}(p) = tempAmpSide(iFreq,i);
%     highLow(p) = rcResultStruct_byCondition{i}.projAvg.errA(IM_index,RC,:,1);
%     highHigh(p) = rcResultStruct_byCondition{i}.projAvg.errA(IM_index,RC,:,2);
%     highSide(p) = rcResultStruct_byCondition{i}.projSidebands.avgA(IM_index, RC);
    p = p + 1;
  end
p = 1;
  for i = MedContrast
    medAmp{iFreq}(p) = tempAmp(iFreq,i);
    medSide{iFreq}(p) = tempAmpSide(iFreq,i);
%     medLow(p) = rcResultStruct_byCondition{i}.projAvg.errA(IM_index,RC,:,1);
%     medHigh(p) = rcResultStruct_byCondition{i}.projAvg.errA(IM_index,RC,:,2);
%     medSide(p) = rcResultStruct_byCondition{i}.projSidebands.avgA(IM_index, RC);
    p = p + 1;
  end
p = 1;
  for i = LowContrast
    lowAmp{iFreq}(p) = tempAmp(iFreq,i);
    lowSide{iFreq}(p) = tempAmpSide(iFreq,i);
%     lowLow(p) = rcResultStruct_byCondition{i}.projAvg.errA(IM_index,RC,:,1);
%     lowHigh(p) = rcResultStruct_byCondition{i}.projAvg.errA(IM_index,RC,:,2);
%     lowSide(p) = rcResultStruct_byCondition{i}.projSidebands.avgA(IM_index, RC);
    p = p + 1;
  end

  side{iFreq} = (highSide{iFreq} + medSide{iFreq} + lowSide{iFreq}) / 3;
% side = [highSide + medSide + lowSide] / 3;

%Amplitudes

temp = flip(highAmp{iFreq});
highAmp{iFreq} = [temp(1:end-1), highAmp{iFreq}];
temp = flip(medAmp{iFreq});
medAmp{iFreq} = [temp(1:end-1), medAmp{iFreq}];
temp = flip(lowAmp{iFreq});
lowAmp{iFreq} = [temp(1:end-1), lowAmp{iFreq}];
temp = flip(side{iFreq});
side{iFreq} = [temp(1:end-1), side{iFreq}];


% %Errors
% temp = flip(highLow);
% highLow = [temp(1:end-1), highLow];
% temp = flip(highHigh);
% highHigh = [temp(1:end-1), highHigh];
% temp = flip(medLow);
% medLow = [temp(1:end-1), medLow];
% temp = flip(highHigh);
% medHigh = [temp(1:end-1), medHigh];
% temp = flip(highLow);
% lowLow = [temp(1:end-1), lowLow];
% temp = flip(highHigh);
% lowHigh = [temp(1:end-1), lowHigh];
% % sidebands
% temp = flip(side);
% side = [temp(1:end-1), side];

colors = brewermap(7, 'BuPu');
noiseColors = brewermap(1,'Greens');

figure;
hold on;
myerrorbar(orientationAll, side{iFreq}, 'Color', noiseColors,'Symbol','o','MarkerSize',10);

% bestfit = fitGaussian(orientationAll, highAmp{iFreq});
% std_coh(iFreq,1) = bestfit.std;
% plot(bestfit.fitX,bestfit.fitY, '-', 'Color', colors(7,:),'LineWidth',3); 
myerrorbar(orientationAll, highAmp{iFreq}, 'Color', colors(7,:),'Symbol','o','MarkerSize',10);%, 'yLow', highLow, 'yHigh', highHigh);

% bestfit = fitGaussian(orientationAll, medAmp{iFreq});
% std_coh(iFreq,2) = bestfit.std;
% plot(bestfit.fitX,bestfit.fitY, '-', 'Color', colors(5,:),'LineWidth',3); 
myerrorbar(orientationAll, medAmp{iFreq}, 'Color', colors(5,:),'Symbol','o','MarkerSize',10);%, 'yLow', medLow, 'yHigh', medHigh);

% bestfit = fitGaussian(orientationAll, lowAmp{iFreq});
% std_coh(iFreq,3) = bestfit.std;
% plot(bestfit.fitX,bestfit.fitY, '-', 'Color', colors(3,:),'LineWidth',3); 
myerrorbar(orientationAll, lowAmp{iFreq}, 'Color', colors(3,:),'Symbol','o','MarkerSize',10);

plot(orientationAll, highAmp{iFreq},'-', 'Color', colors(7,:));
plot(orientationAll, medAmp{iFreq},'-', 'Color', colors(5,:));
plot(orientationAll, lowAmp{iFreq},'-', 'Color', colors(3,:));
plot(orientationAll, side{iFreq}, '-', 'Color', noiseColors);
% 
% myerrorbar(orientationAll, highAmp{iFreq}, 'Color', colors(7,:));%, 'yLow', highLow, 'yHigh', highHigh);
% myerrorbar(orientationAll, medAmp{iFreq}, 'Color', colors(5,:));%, 'yLow', medLow, 'yHigh', medHigh);
% myerrorbar(orientationAll, lowAmp{iFreq}, 'Color', colors(3,:));%, 'yLow', lowLow, 'yHigh', lowHigh);
% myerrorbar(orientationAll, side{iFreq}, 'Color', noiseColors);

xlabel('Orientation (Deg)');
ylabel(sprintf('%s EEG Power',freqLabels{iFreq}));
% legend('High','Med','Low','Noise');

end

fwhm_coh = 2*sqrt(2*log(2)) .* std_coh




highAmpComb = sqrt(highAmp{1}.^2 + highAmp{2}.^2 + highAmp{3}.^2);
medAmpComb = sqrt(medAmp{1}.^2 + medAmp{2}.^2 + medAmp{3}.^2);
lowAmpComb = sqrt(lowAmp{1}.^2 + lowAmp{2}.^2 + lowAmp{3}.^2);
sideComb = sqrt(side{1}.^2 + side{2}.^2 + side{3}.^2);


highAmpComb = sqrt(highAmp{1}.^2 + highAmp{2}.^2 );
medAmpComb = sqrt(medAmp{1}.^2 + medAmp{2}.^2 );
lowAmpComb = sqrt(lowAmp{1}.^2 + lowAmp{2}.^2 );
sideComb = sqrt(side{1}.^2 + side{2}.^2 );


colors = brewermap(7, 'BuPu');
noiseColors = brewermap(1,'Greens');

figure;
hold on;
plot(orientationAll, highAmpComb,'-', 'Color', colors(7,:));
plot(orientationAll, medAmpComb,'-', 'Color', colors(5,:));
plot(orientationAll, lowAmpComb,'-', 'Color', colors(3,:));
plot(orientationAll, sideComb, '-', 'Color', noiseColors);

myerrorbar(orientationAll, highAmpComb, 'Color', colors(7,:));%, 'yLow', highLow, 'yHigh', highHigh);
myerrorbar(orientationAll, medAmpComb, 'Color', colors(5,:));%, 'yLow', medLow, 'yHigh', medHigh);
myerrorbar(orientationAll, lowAmpComb, 'Color', colors(3,:));%, 'yLow', lowLow, 'yHigh', lowHigh);
myerrorbar(orientationAll, sideComb, 'Color', noiseColors);

xlabel('Orientation (Deg)');
ylabel(sprintf('Combined IM EEG Power'));
legend('High','Med','Low','Noise');

