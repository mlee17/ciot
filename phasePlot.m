% subjectList = {'nl-1624','nl-jg','nl-1928','nl-0055','nl-1909','nl-1908','nl-0051','nl-0034','nl-myl','nl-rta','nl-0057','nl-2274','nl-2275'};
% 
% dataDir = '~/data/orientationtuning';
% goodSubjList = {'NL0034_20190726_1755', 'rta_20210730_1409', 'myl_20210730_1548',...
%     'jg_20180309_1321', 'nl-1909_20190801_1508', 'nl-1908_20190816_1730', 'nl-0057_20211011_1016', 'nl-2274_20211012_1155','nl-2275_20211013_1526'};
% badSubjList = {'NL0055_20190809_1139', 'nl-1624_20180313_1554', 'nl-1928_20190806_1602', 'NL0051_20190820_1121'};
% matfileDir = 'Exp_MATL_HCN_128_Avg';

subjectList = {'nl-1624','nl-jg','nl-1928','nl-0055','nl-1909','nl-1908','nl-0051','nl-0034','nl-myl','nl-rta','nl-0057','nl-2274','nl-2275',...
    'nl-2277','nl-2278','nl-2279','nl-2280','nl-2281','nl-2282', 'nl-2283','nl-2285', 'nl-2287', 'nl-2300','nl-1539','nl-2113','nl-2126','nl-2215'};
fileDir = '~/data/orientationRCA/snr';
dataDir = '~/data/orientationtuning';
matfileDir = 'Exp_MATL_HCN_128_Avg';

subjectList = {'nl-1624','nl-jg','nl-1928','nl-0055','nl-1909','nl-1908','nl-0051','nl-0034','nl-myl','nl-rta','nl-0057','nl-2274','nl-2275',...
    'nl-2277','nl-2278','nl-2279','nl-2280','nl-2281','nl-2282', 'nl-2283','nl-2285', 'nl-2287', 'nl-2300','nl-1539','nl-2113','nl-2126','nl-2215'};
subjectList_copy = {'nl-1624','nl-jg','nl-1928','nl-0055','nl-1909','nl-1908','nl-0051','nl-0034','nl-myl','nl-rta','nl-0057','nl-2274','nl-2275',...
    'nl-2277','nl-2278','nl-2279','nl-2280','nl-2281','nl-2282', 'nl-2283','nl-2285', 'nl-2287', 'nl-2300','nl-1539','nl-2113','nl-2126','nl-2215'};

dataDir = '~/data/orientationtuning';
matfileDir = 'Exp_MATL_HCN_128_Avg';

subjectList = subjectList(rankInd);
subjectList = subjectList(isOver2);

nSubj= length(subjectList);

target_frequencies = [8,2];%[8, 2, 4]; %1f1+1f2 -1f1+1f2 -2f1+2f2
 
for s = 1:length(subjectList)
    thisFolder = dir(fullfile(dataDir, [subjectList{s},'*']));
    thisDir = fullfile(dataDir, thisFolder.name, matfileDir);
        
    [output_cos{s}, output_sin{s}, output_amp{s}, output_side_cos{s}, output_side_sin{s}] = getPhases(thisDir, target_frequencies);
    
end
% for s = 1:length(goodSubjList)
%     thisDir = fullfile(dataDir, goodSubjList{s}, matfileDir);
%     
%     [output_cos{s}, output_sin{s}, output_amp{s}, output_side_cos{s}, output_side_sin{s}] = getPhases(thisDir, target_frequencies);
%     
% end
% 
% for s = 1:length(badSubjList)
%     thisDir = fullfile(dataDir, badSubjList{s}, matfileDir);
%     
%     [output_cos2{s}, output_sin2{s}, output_amp2{s}, output_side_cos2{s}, output_side_sin2{s}] = getPhases(thisDir, target_frequencies);
%     
% end
% 
% goodColors = brewermap(13,'*BuGn');
% goodColors = goodColors(3:8,:);
% badColors = brewermap(4,'OrRd');
% 
% freqLabels = {'1F1+1F2', '-1F1+1F2', '-2F1+2F2'};
% 
% figure;
% for iFreq = 1:length(target_frequencies)
%     subplot(1,3,iFreq);
%     hold on;
%     for s = 1:length(goodSubjList)
%         plot([0 output_cos{s}(iFreq,1)], [0 output_sin{s}(iFreq,1)], 'color', goodColors(s,:));
%         hold on;
%     end
%     for s = 1:length(badSubjList)
%         plot([0 output_cos2{s}(iFreq,1)], [0 output_sin2{s}(iFreq,1)], 'color', badColors(s,:));
%         hold on;
%     end
%     xlabel('cos');
%     ylabel('sin');
%     yaxis(-1,1);
%     xaxis(-1,0.5);
%     hline(0);
%     vline(0);
%     title(freqLabels{iFreq});
%     set(gca,'FontSize',13);
% end

%%%%%%
%% Coherent Averaging

tempCos = zeros(2,27); %nFreq, nCond
tempSin = zeros(2,27);
sideCos = zeros(2,27);
sideSin = zeros(2,27);
for s = 1:20%length(subjectList)
    tempCos = tempCos + output_cos{s};
    tempSin = tempSin + output_sin{s};
    sideCos = sideCos + output_side_cos{s};
    sideSin = sideSin + output_side_sin{s};
end

% for s = 1:length(goodSubjList)
%     tempCos = tempCos + output_cos{s};
%     tempSin = tempSin + output_sin{s};
%     sideCos = sideCos + output_side_cos{s};
%     sideSin = sideSin + output_side_sin{s};
% end
% 
% for s = 1:length(badSubjList)
%     tempCos = tempCos + output_cos2{s};
%     tempSin = tempSin + output_sin2{s};
%     sideCos = sideCos + output_side_cos2{s};
%     sideSin = sideSin + output_side_sin2{s};
% end

% for s = 1:length(goodSubjList)
%     tempCos = tempCos + squeeze(mean(output_cos{s},2));
%     tempSin = tempSin + squeeze(mean(output_sin{s},2));
%     sideCos = sideCos + squeeze(mean(output_side_cos{s},2));
%     sideSin = sideSin + squeeze(mean(output_side_sin{s},2));
% end

% tempCos = tempCos / length(goodSubjList);
% tempSin = tempSin / length(goodSubjList);
% sideCos = sideCos / length(goodSubjList);
% sideSin = sideSin / length(goodSubjList);

% tempCos = tempCos / (length(goodSubjList)+length(badSubjList));
% tempSin = tempSin / (length(goodSubjList)+length(badSubjList));
% sideCos = sideCos / (length(goodSubjList)+length(badSubjList));
% sideSin = sideSin / (length(goodSubjList)+length(badSubjList));

tempCos = tempCos / 20;%(length(subjectList));
tempSin = tempSin / 20;%(length(subjectList));
sideCos = sideCos / 20;%(length(subjectList));
sideSin = sideSin / 20;%(length(subjectList));


% for s = 1:length(goodSubjList)
%     tempCos = tempCos + output_cos{s};
%     tempSin = tempSin + output_sin{s};
%     
%     sideCos = sideCos + output_side_cos{s};
%     sideSin = sideSin + output_side_sin{s};
%     
% end
% for s = 1:length(badSubjList)
%     tempCos = tempCos + output_cos2{s};
%     tempSin = tempSin + output_sin2{s};
%     
%     sideCos = sideCos + output_side_cos2{s};
%     sideSin = sideSin + output_side_sin2{s};
% end
% tempCos = tempCos / (length(goodSubjList)+length(badSubjList));
% tempSin = tempSin / (length(goodSubjList)+length(badSubjList));
% sideCos = sideCos / (length(goodSubjList)+length(badSubjList));
% sideSin = sideSin / (length(goodSubjList)+length(badSubjList));

tempAmp = sqrt(tempCos.^2 + tempSin.^2);
tempAmpSide = sqrt(sideCos.^2 + sideSin.^2);

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

bestfit = fitGaussian(orientationAll, highAmp{iFreq});
std_coh(iFreq,1) = bestfit.std;
plot(bestfit.fitX,bestfit.fitY, '-', 'Color', colors(7,:),'LineWidth',3); 
myerrorbar(orientationAll, highAmp{iFreq}, 'Color', colors(7,:),'Symbol','o','MarkerSize',10);%, 'yLow', highLow, 'yHigh', highHigh);

bestfit = fitGaussian(orientationAll, medAmp{iFreq});
std_coh(iFreq,2) = bestfit.std;
plot(bestfit.fitX,bestfit.fitY, '-', 'Color', colors(5,:),'LineWidth',3); 
myerrorbar(orientationAll, medAmp{iFreq}, 'Color', colors(5,:),'Symbol','o','MarkerSize',10);%, 'yLow', medLow, 'yHigh', medHigh);

bestfit = fitGaussian(orientationAll, lowAmp{iFreq});
std_coh(iFreq,3) = bestfit.std;
plot(bestfit.fitX,bestfit.fitY, '-', 'Color', colors(3,:),'LineWidth',3); 
myerrorbar(orientationAll, lowAmp{iFreq}, 'Color', colors(3,:),'Symbol','o','MarkerSize',10);

% plot(orientationAll, highAmp{iFreq},'-', 'Color', colors(7,:));
% plot(orientationAll, medAmp{iFreq},'-', 'Color', colors(5,:));
% plot(orientationAll, lowAmp{iFreq},'-', 'Color', colors(3,:));
% plot(orientationAll, side{iFreq}, '-', 'Color', noiseColors);

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





