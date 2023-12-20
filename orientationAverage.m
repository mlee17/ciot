% subjectList = {'nl-1624','nl-jg','nl-1928','nl-0055','nl-1909','nl-1908','nl-0051','nl-0034','nl-myl','nl-rta'};
% 
% goodSubjList = {'NL0034_20190726_1755', 'rta_20210730_1409', 'myl_20210730_1548',...
%     'jg_20180309_1321', 'nl-1909_20190801_1508', 'nl-1908_20190816_1730'};
% badSubjList = {'NL0055_20190809_1139', 'nl-1624_20180313_1554', 'nl-1928_20190806_1602', 'NL0051_20190820_1121'};
subjectList = {'nl-1624','nl-jg','nl-1928','nl-0055','nl-1909','nl-1908','nl-0051','nl-0034','nl-myl','nl-rta','nl-0057','nl-2274','nl-2275',...
    'nl-2277','nl-2278','nl-2279','nl-2280','nl-2281','nl-2282', 'nl-2283','nl-2285', 'nl-2287', 'nl-2300','nl-1539','nl-2113','nl-2126','nl-2215',...
    'nl-2329','nl-2331','nl-2332','nl-2335','nl-2336','nl-2343','nl-2344','nl-2346','nl-2347','nl-2349','nl-2350','nl-2351','nl-2353', 'nl-2374',...
    'nl-2378', 'nl-2379','nl-2381', 'nl-2383'};
subjectList_copy = {'nl-1624','nl-jg','nl-1928','nl-0055','nl-1909','nl-1908','nl-0051','nl-0034','nl-myl','nl-rta','nl-0057','nl-2274','nl-2275',...
    'nl-2277','nl-2278','nl-2279','nl-2280','nl-2281','nl-2282', 'nl-2283','nl-2285', 'nl-2287', 'nl-2300','nl-1539','nl-2113','nl-2126','nl-2215',...
    'nl-2329','nl-2331','nl-2332','nl-2335','nl-2336','nl-2343','nl-2344','nl-2346','nl-2347','nl-2349','nl-2350','nl-2351','nl-2353', 'nl-2374',...
    'nl-2378', 'nl-2379','nl-2381', 'nl-2383'};

subjectList = {'nl-1624','nl-jg','nl-1928','nl-0055','nl-1909','nl-1908','nl-0051','nl-0034','nl-myl','nl-rta','nl-0057','nl-2274','nl-2275',...
    'nl-2277','nl-2278','nl-2279','nl-2280','nl-2281','nl-2282', 'nl-2283','nl-2285', 'nl-2287', 'nl-2300','nl-1539','nl-2113','nl-2126','nl-2215',...
    'nl-2329','nl-2331','nl-2332','nl-2286'};%,'nl-2335','nl-2336','nl-2343','nl-2344','nl-2346','nl-2347','nl-2349','nl-2350','nl-2351','nl-2353', 'nl-2374',...
    %'nl-2378', 'nl-2379','nl-2381', 'nl-2383'};
    
% task = [zeroes(0,27),0,0,0,ones(0,10)];
% taskList = [zeros(1,27),0,0,0,ones(1,15)];
% taskList_copy = [zeros(1,27),0,0,0,ones(1,15)];
taskList = [zeros(1,31)];

dataDir = '~/data/orientationtuning';
matfileDir = 'Exp_MATL_HCN_128_Avg';

channel_to_visualize = 75;
target_freq = 4;
freq_label = 'F1+F2';

[orderedSNR, rankInd] = sort(meanSNR(1:nSubj,5), 'descend');
isOver2 = orderedSNR>2;

subjectList = subjectList(rankInd);
subjectList = subjectList(isOver2);
taskList = task(rankInd);
taskList = taskList(isOver2);

nSubj= length(subjectList);

% subjectList = {'nl-2336','nl-2350','nl-2346','nl-2374','nl-2383'};%,'nl-2332'};
% taskList = [1,1,1,1,1];
% nSubj= length(subjectList);
for s = 1:nSubj
    thisFolder = dir(fullfile(dataDir, [subjectList{s},'*']));
    thisDir = fullfile(dataDir, thisFolder.name, matfileDir);
    subj_label = subjectList{s};
    
    [high_contrasts(s,:), medium_contrasts(s,:), low_contrasts(s,:), average_contrasts(s,:)] = ...
        orientationPlotSubj(thisDir, channel_to_visualize, target_freq, freq_label, subj_label);
end

for s = 1:nSubj
    thisFolder = dir(fullfile(dataDir, [subjectList{s},'*']));
    thisDir = fullfile(dataDir, thisFolder.name, matfileDir);
    subj_label = subjectList{s};
    
    [high_combined(s,:), medium_combined(s,:), low_combined(s,:), average_combined(s,:)] = ...
        combinedIMSubj(thisDir, subj_label);
end

for s = 1:nSubj
    thisFolder = dir(fullfile(dataDir, [subjectList{s},'*']));
    thisDir = fullfile(dataDir, thisFolder.name, matfileDir);
    subj_label = subjectList{s};
    
    [high_contrastsOcc(s,:), medium_contrastsOcc(s,:), low_contrastsOcc(s,:), average_contrastsOcc(s,:), mean_singleOcc(:,:,s), mean_single_averageOcc(:,s)] = ...
        occipitalPlotSubj(thisDir, subj_label, taskList(s));
end


for s = 1:nSubj
    thisFolder = dir(fullfile(dataDir, [subjectList{s},'*']));
    thisDir = fullfile(dataDir, thisFolder.name, matfileDir);
    subj_label = subjectList{s};
    
    [high_contrastsOcc3(s,:), medium_contrastsOcc3(s,:), low_contrastsOcc3(s,:), average_contrastsOcc3(s,:), mean_singleOcc3(:,:,s), mean_single_averageOcc3(:,s)] = ...
        occipitalPlotSubj(thisDir, subj_label, taskList(s));
end

for s = 1:nSubj
    thisFolder = dir(fullfile(dataDir, [subjectList{s},'*']));
    thisDir = fullfile(dataDir, thisFolder.name, matfileDir);
    subj_label = subjectList{s};
    
    [high_contrastsOcc2(s,:), medium_contrastsOcc2(s,:), low_contrastsOcc2(s,:), average_contrastsOcc2(s,:), mean_singleOcc2(:,:,s), mean_single_averageOcc2(:,s)] = ...
        occipitalPlot_diff(thisDir, subj_label, taskList(s));
end


for s = 1:10
[high_contrasts(s,:), medium_contrasts(s,:), low_contrasts(s,:), average_contrasts(s,:)] = combinedIM;
end

for s = 1:10
[high_contrasts(s,:), medium_contrasts(s,:), low_contrasts(s,:), average_contrasts(s,:)] = orientationPlot(2, '1-27', 75, {'1-7','10-16','19-25'}, 8, 'myOrientationPlot');
end

for s = 1:10
[high_contrasts(s,:), medium_contrasts(s,:), low_contrasts(s,:), average_contrasts(s,:)] = occipitalPlot;
end

for s = 1:10
    [high_contrasts(s,:), medium_contrasts(s,:), low_contrasts(s,:)] = significantSNRPlot(2, '1-27', {'1-7','10-16','19-25'}, 8);
end

for s = 1:10
[high_contrasts(s,:), medium_contrasts(s,:), low_contrasts(s,:), average_contrasts(s,:)] = orientationPlot(2, '1-27', 76, {'1-7','10-16','19-25'}, 2, 'myOrientationPlot');
end

for s = 1:10
[high_contrasts(s,:), medium_contrasts(s,:), low_contrasts(s,:), average_contrasts(s,:)] = orientationPlot(2, '1-27', 76, {'1-7','10-16','19-25'}, 16, 'myOrientationPlot');
end

orientationOffsets = [0 7 15 30 45 75 90];
orientationAll = [-90 -75 -45 -30 -15 -7 0 7 15 30 45 75 90];

highAve1 = mean(high_contrasts, 1);
medAve1 = mean(medium_contrasts, 1);
lowAve1 = mean(low_contrasts, 1);
highErr1 = std(high_contrasts, 0, 1) / sqrt(size(high_contrasts,1));
medErr1 = std( medium_contrasts, 0, 1) / sqrt(size(medium_contrasts,1));
lowErr1 = std(low_contrasts, 0, 1) / sqrt(size(low_contrasts,1));

noiseAve1 = mean(average_contrasts, 1);
noiseErr1 = std(average_contrasts, 0, 1) / sqrt(size(average_contrasts,1));

highAve2 = mean(high_contrasts, 1);
medAve2 = mean(medium_contrasts, 1);
lowAve2 = mean(low_contrasts, 1);
highErr2 = std(high_contrasts, 0, 1) / sqrt(size(high_contrasts,1));
medErr2 = std( medium_contrasts, 0, 1) / sqrt(size(medium_contrasts,1));
lowErr2 = std(low_contrasts, 0, 1) / sqrt(size(low_contrasts,1));

noiseAve2 = mean(average_contrasts, 1);
noiseErr2 = std(average_contrasts, 0, 1) / sqrt(size(average_contrasts,1));

%%%%%
highAveOcc = mean(high_contrastsOcc, 1);
medAveOcc = mean(medium_contrastsOcc, 1);
lowAveOcc = mean(low_contrastsOcc, 1);
highErrOcc = std(high_contrastsOcc, 0, 1) / sqrt(size(high_contrastsOcc,1));
medErrOcc = std( medium_contrastsOcc, 0, 1) / sqrt(size(medium_contrastsOcc,1));
lowErrOcc = std(low_contrastsOcc, 0, 1) / sqrt(size(low_contrastsOcc,1));
lowErrOcc_std = std(low_contrastsOcc, 0, 1);

noiseAveOcc = mean(average_contrastsOcc, 1);
noiseErrOcc = std(average_contrastsOcc, 0, 1) / sqrt(size(average_contrastsOcc,1));

singleAveOcc = mean(mean_singleOcc,3);
singleErrOcc = std(mean_singleOcc, 0, 3) / sqrt(size(mean_singleOcc, 3));
singleNoiseOcc = mean(mean_single_averageOcc,2);
singleNoiseErrOcc = std(mean_single_averageOcc, 0, 2) / sqrt(size(mean_single_averageOcc,2));

%%%%%
highAveOcc2 = mean(high_contrastsOcc2, 1);
medAveOcc2 = mean(medium_contrastsOcc2, 1);
lowAveOcc2 = mean(low_contrastsOcc2, 1);
highErrOcc2 = std(high_contrastsOcc2, 0, 1) / sqrt(size(high_contrastsOcc2,1));
medErrOcc2 = std( medium_contrastsOcc2, 0, 1) / sqrt(size(medium_contrastsOcc2,1));
lowErrOcc2 = std(low_contrastsOcc2, 0, 1) / sqrt(size(low_contrastsOcc2,1));
lowErrOcc2_std = std(low_contrastsOcc2, 0, 1);

noiseAveOcc2 = mean(average_contrastsOcc2, 1);
noiseErrOcc2 = std(average_contrastsOcc2, 0, 1) / sqrt(size(average_contrastsOcc2,1));

singleAveOcc2 = mean(mean_singleOcc2,3);
singleErrOcc2 = std(mean_singleOcc2, 0, 3) / sqrt(size(mean_singleOcc2, 3));
singleNoiseOcc2 = mean(mean_single_averageOcc2,2);
singleNoiseErrOcc2 = std(mean_single_averageOcc2, 0, 2) / sqrt(size(mean_single_averageOcc2,2));


% 
% highAve2 = mean(high_contrasts(7:10,:), 1);
% medAve2 = mean(medium_contrasts(7:10,:), 1);
% lowAve2 = mean(low_contrasts(7:10,:), 1);
% highErr2 = std(high_contrasts(7:10,:), 0, 1) / sqrt(4);%(size(high_contrasts,1)-2);
% medErr2 = std( medium_contrasts(7:10,:), 0, 1) / sqrt(4);%(size(medium_contrasts,1)-2);
% lowErr2 = std(low_contrasts(7:10,:), 0, 1) / sqrt(4);%(size(low_contrasts,1)-2);
% 
% noiseAve2 = mean(average_contrasts(7:10,:), 1);
% noiseErr2 = std(average_contrasts(7:10,:), 0, 1) / sqrt(4);%(size(average_contrasts,1));



colors = brewermap(7, 'BuPu');
noiseColors = brewermap(1,'Greens');

figure;
hold on;

plot(orientationAll, highAve1,'-', 'Color', colors(7,:));
plot(orientationAll, medAve1,'-', 'Color', colors(4,:));
plot(orientationAll, lowAve1,'-', 'Color', colors(3,:));
plot(orientationAll, noiseAve1, '-', 'Color', noiseColors);

legend('High','Med','Low','Noise');

myerrorbar(orientationAll, highAve1, 'Color', colors(7,:), 'yError', highErr1,'Symbol','o','MarkerSize',10);
myerrorbar(orientationAll, medAve1, 'Color', colors(4,:), 'yError', medErr1,'Symbol','o','MarkerSize',10);
myerrorbar(orientationAll, lowAve1, 'Color', colors(3,:), 'yError', lowErr1,'Symbol','o','MarkerSize',10);
myerrorbar(orientationAll, noiseAve1, 'Color', noiseColors, 'yError', noiseErr1,'Symbol','o','MarkerSize',10);

xlabel('Orientation (Deg)');
ylabel('-F1+F2 EEG Power');
ylabel('Combined IM EEG Power');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
hold on;
plot(orientationAll, highAve2,'-', 'Color', colors(7,:));
plot(orientationAll, medAve2,'-', 'Color', colors(5,:));
plot(orientationAll, lowAve2,'-', 'Color', colors(3,:));
plot(orientationAll, noiseAve2, '-', 'Color', noiseColors);

legend('High','Med','Low','Noise');

myerrorbar(orientationAll, highAve2, 'Color', colors(7,:), 'yError', highErr2,'Symbol','o','MarkerSize',10);
myerrorbar(orientationAll, medAve2, 'Color', colors(5,:), 'yError', medErr2,'Symbol','o','MarkerSize',10);
myerrorbar(orientationAll, lowAve2, 'Color', colors(3,:), 'yError', lowErr2,'Symbol','o','MarkerSize',10);
myerrorbar(orientationAll, noiseAve2, 'Color', noiseColors, 'yError', noiseErr2,'Symbol','o','MarkerSize',10);


xlabel('Orientation (Deg)');
ylabel('F1+F2 EEG Power');
ylabel('Combined IM EEG Power');

%%%%%%%%%%%%%%%%%%

colors = brewermap(7, 'BuPu');
noiseColors = brewermap(1,'Greens');

figure;
hold on;

plot(orientationAll, noiseAveOcc, '-', 'Color', noiseColors);
plot(orientationAll, highAveOcc,'-', 'Color', colors(7,:));
plot(orientationAll, medAveOcc,'-', 'Color', colors(5,:));
plot(orientationAll, lowAveOcc,'-', 'Color', colors(3,:));

legend('Noise','High','Med','Low');

myerrorbar(orientationAll, noiseAveOcc, 'Color', noiseColors, 'yError', noiseErrOcc,'Symbol','o','MarkerSize',10);
myerrorbar(orientationAll, highAveOcc, 'Color', colors(7,:), 'yError', highErrOcc,'Symbol','o','MarkerSize',10);
myerrorbar(orientationAll, medAveOcc, 'Color', colors(5,:), 'yError', medErrOcc,'Symbol','o','MarkerSize',10);
myerrorbar(orientationAll, lowAveOcc, 'Color', colors(3,:), 'yError', lowErrOcc ,'Symbol','o','MarkerSize',10);

% myerrorbar(110, singleAveOcc(1,1)', 'Color', colors(7,:), 'yError', singleErrOcc(1,1),'Symbol','o', 'MarkerSize',10);
% myerrorbar(110, singleAveOcc(1,2), 'Color', colors(5,:), 'yError', singleErrOcc(1,2),'Symbol','o', 'MarkerSize',10);
% myerrorbar(110, singleAveOcc(1,3), 'Color', colors(3,:), 'yError', singleErrOcc(1,3),'Symbol','o', 'MarkerSize',10);
% myerrorbar(130, singleAveOcc(2,1), 'Color', colors(7,:), 'yError', singleErrOcc(2,1),'Symbol','o', 'MarkerSize',10);
% myerrorbar(130, singleAveOcc(2,2), 'Color', colors(5,:), 'yError', singleErrOcc(2,2),'Symbol','o', 'MarkerSize',10);
% myerrorbar(130, singleAveOcc(2,3), 'Color', colors(3,:), 'yError', singleErrOcc(2,3),'Symbol','o', 'MarkerSize',10);
%     
% myerrorbar(110, singleNoiseOcc(1), 'Color', noiseColors, 'yError',singleNoiseErrOcc(1),'Symbol','o', 'MarkerSize', 10);
% myerrorbar(130, singleNoiseOcc(2), 'Color', noiseColors, 'yError',singleNoiseErrOcc(2),'Symbol','o', 'MarkerSize', 10);

% xlim([-100,130]);

xlabel('Orientation Offset (Deg)');
ylabel('f1+f2 Amplitude (µV)');
% ylabel('Combined IM EEG Power');

figure;
hold on;

plot(orientationAll, noiseAveOcc2, '-', 'Color', noiseColors);
plot(orientationAll, highAveOcc2,'-', 'Color', colors(7,:));
plot(orientationAll, medAveOcc2,'-', 'Color', colors(5,:));
plot(orientationAll, lowAveOcc2,'-', 'Color', colors(3,:));

% legend('High','Med','Low','Noise');

myerrorbar(orientationAll, noiseAveOcc2, 'Color', noiseColors, 'yError', noiseErrOcc2,'Symbol','o','MarkerSize',10);
myerrorbar(orientationAll, highAveOcc2, 'Color', colors(7,:), 'yError', highErrOcc2,'Symbol','o','MarkerSize',10);
myerrorbar(orientationAll, medAveOcc2, 'Color', colors(5,:), 'yError', medErrOcc2,'Symbol','o','MarkerSize',10);
myerrorbar(orientationAll, lowAveOcc2, 'Color', colors(3,:), 'yError', lowErrOcc2,'Symbol','o','MarkerSize',10);

xlabel('Orientation Offset (Deg)');
ylabel('f2-f1 Amplitude (µV)');

%%%

highAve_comb = mean(high_combined, 1);
medAve_comb = mean(medium_combined, 1);
lowAve_comb = mean(low_combined, 1);
highErr_comb = std(high_combined, 0, 1) / sqrt(size(high_combined,1));
medErr_comb = std( medium_combined, 0, 1) / sqrt(size(medium_combined,1));
lowErr_comb = std(low_combined, 0, 1) / sqrt(size(low_combined,1));

noiseAve_comb = mean(average_combined, 1);
noiseErr_comb = std(average_combined, 0, 1) / sqrt(size(average_combined,1));

figure;
hold on;

plot(orientationAll, noiseAve_comb, '-', 'Color', noiseColors);
plot(orientationAll, highAve_comb,'-', 'Color', colors(7,:));
plot(orientationAll, medAve_comb,'-', 'Color', colors(5,:));
plot(orientationAll, lowAve_comb,'-', 'Color', colors(3,:));

% legend('High','Med','Low','Noise');

myerrorbar(orientationAll, noiseAve_comb, 'Color', noiseColors, 'yError', noiseErr_comb,'Symbol','o','MarkerSize',10);
myerrorbar(orientationAll, highAve_comb, 'Color', colors(7,:), 'yError', highErr_comb,'Symbol','o','MarkerSize',10);
myerrorbar(orientationAll, medAve_comb, 'Color', colors(5,:), 'yError', medErr_comb,'Symbol','o','MarkerSize',10);
myerrorbar(orientationAll, lowAve_comb, 'Color', colors(3,:), 'yError', lowErr_comb,'Symbol','o','MarkerSize',10);

xlabel('Orientation Offset (Deg)');
ylabel('Combined IM Amplitude (µV)');


