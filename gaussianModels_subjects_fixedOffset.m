subjectList = {'nl-1624','nl-jg','nl-1928','nl-0055','nl-1909','nl-1908','nl-0051','nl-0034','nl-myl','nl-rta','nl-0057','nl-2274','nl-2275',...
    'nl-2277','nl-2278','nl-2279','nl-2280','nl-2281','nl-2282', 'nl-2283','nl-2285', 'nl-2287', 'nl-2300','nl-1539','nl-2113','nl-2126','nl-2215'};
subjectList_copy = {'nl-1624','nl-jg','nl-1928','nl-0055','nl-1909','nl-1908','nl-0051','nl-0034','nl-myl','nl-rta','nl-0057','nl-2274','nl-2275',...
    'nl-2277','nl-2278','nl-2279','nl-2280','nl-2281','nl-2282', 'nl-2283','nl-2285', 'nl-2287', 'nl-2300','nl-1539','nl-2113','nl-2126','nl-2215'};

dataDir = '~/data/orientationtuning';
matfileDir = 'Exp_MATL_HCN_128_Avg';

channel_to_visualize = 75;
target_freq = 8;
freq_label = 'F1+F2';

% [orderedSNR, rankInd] = sort(meanSNR(1:nSubj,10), 'descend');
% isOver2 = orderedSNR>2;

subjectList = subjectList(rankInd);
subjectList = subjectList(isOver2);

nSubj= length(subjectList);
F_subj_fixed = zeros(nSubj, 1);
orientationAll = [-90 -75 -45 -30 -15 -7 0 7 15 30 45 75 90];
for s = 1:nSubj
    thisFolder = dir(fullfile(dataDir, [subjectList{s},'*']));
    thisDir = fullfile(dataDir, thisFolder.name, matfileDir);
    subj_label = subjectList{s};
    
    [high_contrastsOcc(s,:), medium_contrastsOcc(s,:), low_contrastsOcc(s,:), average_contrastsOcc(s,:), mean_singleOcc(:,:,s), mean_single_averageOcc(:,s)] = ...
        occipitalPlotSubj(thisDir, subj_label);
    
    [bestfit_full, bestfit_amp, F_subj_fixed(s,1), Fr_subj(s,1), AIC_full_subj(s,1), AIC_amp_subj(s,1)] = ...
        compareGaussianModels_fixedOffset(orientationAll, high_contrastsOcc(s,:), medium_contrastsOcc(s,:), low_contrastsOcc(s,:));
    
    modelComparisonPlot(orientationAll, high_contrastsOcc(s,:), medium_contrastsOcc(s,:), low_contrastsOcc(s,:), average_contrastsOcc(s,:),...
        [], [], [], [], bestfit_full, bestfit_amp, F_subj_fixed(s,1), subj_label, 'F1+F2')
    
end
figure;
F_sorted = sort(F_subj_fixed);
plot(F_sorted,ones(length(F_sorted),1)+(rand(length(F_sorted),1)/10), 'o','color',[1 1 1], 'markerFaceColor',[0.2 0.2 0.2],'markerSize',15);
yaxis(0,2);
vline(3.305,'--');
xlabel('F-value');
set(gca,'ytick',[]);
title('f1+f2');

AIC_diff_subj = AIC_amp_subj - AIC_full_subj;
figure;
AIC_sorted = sort(AIC_diff_subj);
plot(AIC_sorted, ones(length(AIC_sorted),1)+(rand(length(F_sorted),1)/10), 'o','color',[1 1 1], 'markerFaceColor',[0.2 0.2 0.2],'markerSize',15);
yaxis(0,2);
vline(0,'--');
xlabel('AIC diff (amp-full)');
set(gca,'ytick',[]);
title('f1+f2');

figure; hold on;
plot(F_subj_fixed, AIC_diff_subj,'o');
vline(3.305,'--');
hline(0,'--');
xlabel('F'); ylabel('AIC diff');


%%%%%
highAveOcc = mean(high_contrastsOcc, 1);
medAveOcc = mean(medium_contrastsOcc, 1);
lowAveOcc = mean(low_contrastsOcc, 1);
highErrOcc = std(high_contrastsOcc, 0, 1) / sqrt(size(high_contrastsOcc,1));
medErrOcc = std( medium_contrastsOcc, 0, 1) / sqrt(size(medium_contrastsOcc,1));
lowErrOcc = std(low_contrastsOcc, 0, 1) / sqrt(size(low_contrastsOcc,1));
noiseAveOcc = mean(average_contrastsOcc, 1);
noiseErrOcc = std(average_contrastsOcc, 0, 1) / sqrt(size(average_contrastsOcc,1));
%%%
% 
% highAveOcc = mean(high_contrastsOcc(1:10,:), 1);
% medAveOcc = mean(medium_contrastsOcc(1:10,:), 1);
% lowAveOcc = mean(low_contrastsOcc(1:10,:), 1);
% highErrOcc = std(high_contrastsOcc(1:10,:), 0, 1) / sqrt(size(high_contrastsOcc(1:10,:),1));
% medErrOcc = std( medium_contrastsOcc(1:10,:), 0, 1) / sqrt(size(medium_contrastsOcc(1:10,:),1));
% lowErrOcc = std(low_contrastsOcc(1:10,:), 0, 1) / sqrt(size(low_contrastsOcc(1:10,:),1));
% noiseAveOcc = mean(average_contrastsOcc(1:10,:), 1);
% noiseErrOcc = std(average_contrastsOcc(1:10,:), 0, 1) / sqrt(size(average_contrastsOcc(1:10,:),1));
% 
% highAveOcc = mean(high_contrastsOcc(11:20,:), 1);
% medAveOcc = mean(medium_contrastsOcc(11:20,:), 1);
% lowAveOcc = mean(low_contrastsOcc(11:20,:), 1);
% highErrOcc = std(high_contrastsOcc(11:20,:), 0, 1) / sqrt(size(high_contrastsOcc(11:20,:),1));
% medErrOcc = std( medium_contrastsOcc(11:20,:), 0, 1) / sqrt(size(medium_contrastsOcc(11:20,:),1));
% lowErrOcc = std(low_contrastsOcc(11:20,:), 0, 1) / sqrt(size(low_contrastsOcc(11:20,:),1));
% noiseAveOcc = mean(average_contrastsOcc(11:20,:), 1);
% noiseErrOcc = std(average_contrastsOcc(11:20,:), 0, 1) / sqrt(size(average_contrastsOcc(11:20,:),1));


[bestfit_full, bestfit_amp, F_obt_fixed, Fr, AIC_full, AIC_amp] =compareGaussianModels_fixedOffset(orientationAll, highAveOcc, medAveOcc, lowAveOcc);
figure; 
subplot(1,2,1); hold on;
myerrorbar(orientationAll, noiseAveOcc, 'Color', noiseColors, 'yError', noiseErrOcc,'Symbol','o','MarkerSize',10);
plot(bestfit_full.fitX, bestfit_full.fitY_high, '-', 'Color', colors(7,:),'LineWidth',3); 
plot(bestfit_full.fitX, bestfit_full.fitY_med, '-', 'Color', colors(5,:),'LineWidth',3); 
plot(bestfit_full.fitX, bestfit_full.fitY_low, '-', 'Color', colors(3,:),'LineWidth',3); 
myerrorbar(orientationAll, highAveOcc, 'Color', colors(7,:), 'yError', highErrOcc,'Symbol','o','MarkerSize',10);
myerrorbar(orientationAll, medAveOcc, 'Color', colors(5,:), 'yError', medErrOcc,'Symbol','o','MarkerSize',10);
myerrorbar(orientationAll, lowAveOcc, 'Color', colors(3,:), 'yError', lowErrOcc ,'Symbol','o','MarkerSize',10);
set(gca,'FontSize', 13);
xlabel('Orientation Offset (Deg)');
ylabel('f1+f2  Amplitude (�V)');
title(sprintf('Full Model Fobt = %0.4f', F_obt_fixed));
subplot(1,2,2); hold on;
myerrorbar(orientationAll, noiseAveOcc, 'Color', noiseColors, 'yError', noiseErrOcc,'Symbol','o','MarkerSize',10);
plot(bestfit_amp.fitX, bestfit_amp.fitY_high, '-', 'Color', colors(7,:),'LineWidth',3); 
plot(bestfit_amp.fitX, bestfit_amp.fitY_med, '-', 'Color', colors(5,:),'LineWidth',3); 
plot(bestfit_amp.fitX, bestfit_amp.fitY_low, '-', 'Color', colors(3,:),'LineWidth',3); 
myerrorbar(orientationAll, highAveOcc, 'Color', colors(7,:), 'yError', highErrOcc,'Symbol','o','MarkerSize',10);
myerrorbar(orientationAll, medAveOcc, 'Color', colors(5,:), 'yError', medErrOcc,'Symbol','o','MarkerSize',10);
myerrorbar(orientationAll, lowAveOcc, 'Color', colors(3,:), 'yError', lowErrOcc ,'Symbol','o','MarkerSize',10);
set(gca,'FontSize', 13);
xlabel('Orientation Offset (Deg)');
ylabel('f1+f2  Amplitude (�V)');
title(sprintf('Amplitude Model Fobt = %0.4f', F_obt_fixed));

F_obt_fixed
% p
% degrees of freedom = (2,28)

%%
% 
for s = 1:nSubj
    thisFolder = dir(fullfile(dataDir, [subjectList{s},'*']));
    thisDir = fullfile(dataDir, thisFolder.name, matfileDir);
    subj_label = subjectList{s};
    
    [high_contrastsOcc2(s,:), medium_contrastsOcc2(s,:), low_contrastsOcc2(s,:), average_contrastsOcc2(s,:), mean_singleOcc2(:,:,s), mean_single_averageOcc2(:,s)] = ...
        occipitalPlot_diff(thisDir, subj_label, taskList(s));
    
    [bestfit_full, bestfit_amp, F_subj_fixed2(s,1), Fr_subj2(s,1), AIC_full_subj2(s,1), AIC_amp_subj2(s,1)] =...
        compareGaussianModels_fixedOffset(orientationAll, high_contrastsOcc2(s,:), medium_contrastsOcc2(s,:), low_contrastsOcc2(s,:));
    
    modelComparisonPlot(orientationAll, high_contrastsOcc2(s,:), medium_contrastsOcc2(s,:), low_contrastsOcc2(s,:), average_contrastsOcc2(s,:),...
        [], [], [], [], bestfit_full, bestfit_amp, F_subj_fixed2(s,1), subj_label,'F2-F1')
end
figure;
F_sorted = sort(F_subj_fixed2);
plot(F_sorted,ones(length(F_sorted),1)+(rand(length(F_sorted),1)/10), 'o','color',[1 1 1], 'markerFaceColor',[0.2 0.2 0.2],'markerSize',15);
yaxis(0,2);
vline(3.305,'--');
xlabel('F-value');
set(gca,'ytick',[]);
title('f2-f1');

AIC_diff_subj2 = AIC_amp_subj2 - AIC_full_subj2;
figure;
AIC_sorted2 = sort(AIC_diff_subj2);
plot(AIC_sorted2, ones(length(AIC_sorted2),1)+(rand(length(AIC_sorted2),1)/10), 'o','color',[1 1 1], 'markerFaceColor',[0.2 0.2 0.2],'markerSize',15);
yaxis(0,2);
vline(0,'--');
xlabel('AIC diff (amp-full)');
set(gca,'ytick',[]);
title('f2-f1');

figure; hold on;
plot(F_subj_fixed2, AIC_diff_subj2,'o');
vline(3.305,'--');
hline(0,'--');
xlabel('F'); ylabel('AIC diff');

%%%%%
highAveOcc2 = mean(high_contrastsOcc2, 1);
medAveOcc2 = mean(medium_contrastsOcc2, 1);
lowAveOcc2 = mean(low_contrastsOcc2, 1);
highErrOcc2 = std(high_contrastsOcc2, 0, 1) / sqrt(size(high_contrastsOcc2,1));
medErrOcc2 = std( medium_contrastsOcc2, 0, 1) / sqrt(size(medium_contrastsOcc2,1));
lowErrOcc2 = std(low_contrastsOcc2, 0, 1) / sqrt(size(low_contrastsOcc2,1));
noiseAveOcc2 = mean(average_contrastsOcc2, 1);
noiseErrOcc2 = std(average_contrastsOcc2, 0, 1) / sqrt(size(average_contrastsOcc2,1));

%%%
% highAveOcc2 = mean(high_contrastsOcc2(1:10,:), 1);
% medAveOcc2 = mean(medium_contrastsOcc2(1:10,:), 1);
% lowAveOcc2 = mean(low_contrastsOcc2(1:10,:), 1);
% highErrOcc2 = std(high_contrastsOcc2(1:10,:), 0, 1) / sqrt(size(high_contrastsOcc2(1:10,:),1));
% medErrOcc2 = std( medium_contrastsOcc2(1:10,:), 0, 1) / sqrt(size(medium_contrastsOcc2(1:10,:),1));
% lowErrOcc2 = std(low_contrastsOcc2(1:10,:), 0, 1) / sqrt(size(low_contrastsOcc2(1:10,:),1));
% noiseAveOcc2 = mean(average_contrastsOcc2(1:10,:), 1);
% noiseErrOcc2 = std(average_contrastsOcc2(1:10,:), 0, 1) / sqrt(size(average_contrastsOcc2(1:10,:),1));
% 
% highAveOcc2 = mean(high_contrastsOcc2(11:20,:), 1);
% medAveOcc2 = mean(medium_contrastsOcc2(11:20,:), 1);
% lowAveOcc2 = mean(low_contrastsOcc2(11:20,:), 1);
% highErrOcc2 = std(high_contrastsOcc2(11:20,:), 0, 1) / sqrt(size(high_contrastsOcc2(11:20,:),1));
% medErrOcc2 = std( medium_contrastsOcc2(11:20,:), 0, 1) / sqrt(size(medium_contrastsOcc2(11:20,:),1));
% lowErrOcc2 = std(low_contrastsOcc2(11:20,:), 0, 1) / sqrt(size(low_contrastsOcc2(11:20,:),1));
% noiseAveOcc2 = mean(average_contrastsOcc2(11:20,:), 1);
% noiseErrOcc2 = std(average_contrastsOcc2(11:20,:), 0, 1) / sqrt(size(average_contrastsOcc2(11:20,:),1));


[bestfit_full, bestfit_amp, F_obt_fixed2, Fr2, AIC_full2, AIC_amp2] =compareGaussianModels_fixedOffset(orientationAll, highAveOcc2, medAveOcc2, lowAveOcc2);
figure; 
subplot(1,2,1); hold on;
myerrorbar(orientationAll, noiseAveOcc2, 'Color', noiseColors, 'yError', noiseErrOcc2,'Symbol','o','MarkerSize',10);
plot(bestfit_full.fitX, bestfit_full.fitY_high, '-', 'Color', colors(7,:),'LineWidth',3); 
plot(bestfit_full.fitX, bestfit_full.fitY_med, '-', 'Color', colors(5,:),'LineWidth',3); 
plot(bestfit_full.fitX, bestfit_full.fitY_low, '-', 'Color', colors(3,:),'LineWidth',3); 
myerrorbar(orientationAll, highAveOcc2, 'Color', colors(7,:), 'yError', highErrOcc2,'Symbol','o','MarkerSize',10);
myerrorbar(orientationAll, medAveOcc2, 'Color', colors(5,:), 'yError', medErrOcc2,'Symbol','o','MarkerSize',10);
myerrorbar(orientationAll, lowAveOcc2, 'Color', colors(3,:), 'yError', lowErrOcc2 ,'Symbol','o','MarkerSize',10);
set(gca,'FontSize', 13);
xlabel('Orientation Offset (Deg)');
ylabel('f2-f1 Amplitude (�V)');
title(sprintf('Full Model Fobt = %0.4f', F_obt_fixed2));
subplot(1,2,2); hold on;
myerrorbar(orientationAll, noiseAveOcc2, 'Color', noiseColors, 'yError', noiseErrOcc2,'Symbol','o','MarkerSize',10);
plot(bestfit_amp.fitX, bestfit_amp.fitY_high, '-', 'Color', colors(7,:),'LineWidth',3); 
plot(bestfit_amp.fitX, bestfit_amp.fitY_med, '-', 'Color', colors(5,:),'LineWidth',3); 
plot(bestfit_amp.fitX, bestfit_amp.fitY_low, '-', 'Color', colors(3,:),'LineWidth',3); 
myerrorbar(orientationAll, highAveOcc2, 'Color', colors(7,:), 'yError', highErrOcc2,'Symbol','o','MarkerSize',10);
myerrorbar(orientationAll, medAveOcc2, 'Color', colors(5,:), 'yError', medErrOcc2,'Symbol','o','MarkerSize',10);
myerrorbar(orientationAll, lowAveOcc2, 'Color', colors(3,:), 'yError', lowErrOcc2 ,'Symbol','o','MarkerSize',10);
set(gca,'FontSize', 13);
xlabel('Orientation Offset (Deg)');
ylabel('f2-f1 Amplitude (�V)');
title(sprintf('Amplitude Model Fobt = %0.4f', F_obt_fixed2));

F_obt_fixed2
% p2 

%% combined
for s = 1:nSubj
    thisFolder = dir(fullfile(dataDir, [subjectList{s},'*']));
    thisDir = fullfile(dataDir, thisFolder.name, matfileDir);
    subj_label = subjectList{s};
    
    [high_combined(s,:), medium_combined(s,:), low_combined(s,:), average_combined(s,:)] = ...
        combinedIMSubj(thisDir, subj_label, taskList(s));
    
    [bestfit_full, bestfit_amp, F_subj_fixed_combined(s,1)] = compareGaussianModels_fixedOffset(orientationAll, high_combined(s,:), medium_combined(s,:), low_combined(s,:));
    
end
highAve_comb = mean(high_combined, 1);
medAve_comb = mean(medium_combined, 1);
lowAve_comb = mean(low_combined, 1);
highErr_comb = std(high_combined, 0, 1) / sqrt(size(high_combined,1));
medErr_comb = std( medium_combined, 0, 1) / sqrt(size(medium_combined,1));
lowErr_comb = std(low_combined, 0, 1) / sqrt(size(low_combined,1));

noiseAve_comb = mean(average_combined, 1);
noiseErr_comb = std(average_combined, 0, 1) / sqrt(size(average_combined,1));

[bestfit_full, bestfit_amp, F_obt_fixed_combined] = compareGaussianModels_fixedOffset(orientationAll, highAve_comb, medAve_comb, lowAve_comb);
figure;
subplot(1,2,1); hold on;
myerrorbar(orientationAll, noiseAve_comb, 'Color', noiseColors, 'yError', noiseErr_comb,'Symbol','o','MarkerSize',10);
plot(bestfit_full.fitX, bestfit_full.fitY_high, '-', 'Color', colors(7,:),'LineWidth',3); 
plot(bestfit_full.fitX, bestfit_full.fitY_med, '-', 'Color', colors(5,:),'LineWidth',3); 
plot(bestfit_full.fitX, bestfit_full.fitY_low, '-', 'Color', colors(3,:),'LineWidth',3); 
myerrorbar(orientationAll, highAve_comb, 'Color', colors(7,:), 'yError', highErr_comb,'Symbol','o','MarkerSize',10);
myerrorbar(orientationAll, medAve_comb, 'Color', colors(5,:), 'yError', medErr_comb,'Symbol','o','MarkerSize',10);
myerrorbar(orientationAll, lowAve_comb, 'Color', colors(3,:), 'yError', lowErr_comb ,'Symbol','o','MarkerSize',10);
set(gca,'FontSize', 13);
xlabel('Orientation Offset (Deg)');
ylabel('Combined IM Amplitude (�V)');
title(sprintf('Full Model Fobt = %0.4f', F_obt_fixed_combined));
subplot(1,2,2); hold on;
myerrorbar(orientationAll, noiseAve_comb, 'Color', noiseColors, 'yError', noiseErr_comb,'Symbol','o','MarkerSize',10);
plot(bestfit_amp.fitX, bestfit_amp.fitY_high, '-', 'Color', colors(7,:),'LineWidth',3); 
plot(bestfit_amp.fitX, bestfit_amp.fitY_med, '-', 'Color', colors(5,:),'LineWidth',3); 
plot(bestfit_amp.fitX, bestfit_amp.fitY_low, '-', 'Color', colors(3,:),'LineWidth',3); 
myerrorbar(orientationAll, highAve_comb, 'Color', colors(7,:), 'yError', highErr_comb,'Symbol','o','MarkerSize',10);
myerrorbar(orientationAll, medAve_comb, 'Color', colors(5,:), 'yError', medErr_comb,'Symbol','o','MarkerSize',10);
myerrorbar(orientationAll, lowAve_comb, 'Color', colors(3,:), 'yError', lowErr_comb ,'Symbol','o','MarkerSize',10);
set(gca,'FontSize', 13);
xlabel('Orientation Offset (Deg)');
ylabel('Combined IM Amplitude (�V)');
title(sprintf('Amplitude Model Fobt = %0.4f', F_obt_fixed_combined));



%% joint fitting
%%%
[bestfit_full, bestfit_amp, F_obt_fixed_joint] =compareGaussianModels_IM(orientationAll,highAveOcc, medAveOcc, lowAveOcc, highAveOcc2, medAveOcc2, lowAveOcc2);
figure;
subplot(2,2,1); hold on;
myerrorbar(orientationAll, noiseAveOcc, 'Color', noiseColors, 'yError', noiseErrOcc,'Symbol','o','MarkerSize',10);
plot(bestfit_full.fitX, bestfit_full.fitY_high, '-', 'Color', colors(7,:),'LineWidth',3); 
plot(bestfit_full.fitX, bestfit_full.fitY_med, '-', 'Color', colors(5,:),'LineWidth',3); 
plot(bestfit_full.fitX, bestfit_full.fitY_low, '-', 'Color', colors(3,:),'LineWidth',3); 
myerrorbar(orientationAll, highAveOcc, 'Color', colors(7,:), 'yError', highErrOcc,'Symbol','o','MarkerSize',10);
myerrorbar(orientationAll, medAveOcc, 'Color', colors(5,:), 'yError', medErrOcc,'Symbol','o','MarkerSize',10);
myerrorbar(orientationAll, lowAveOcc, 'Color', colors(3,:), 'yError', lowErrOcc ,'Symbol','o','MarkerSize',10);
set(gca,'FontSize', 13);
xlabel('Orientation Offset (Deg)');
ylabel('f1+f2  Amplitude (�V)');
title(sprintf('Full Model Fobt = %0.4f', F_obt_fixed_joint));

subplot(2,2,2); hold on;
myerrorbar(orientationAll, noiseAveOcc, 'Color', noiseColors, 'yError', noiseErrOcc,'Symbol','o','MarkerSize',10);
plot(bestfit_amp.fitX, bestfit_amp.fitY_high, '-', 'Color', colors(7,:),'LineWidth',3); 
plot(bestfit_amp.fitX, bestfit_amp.fitY_med, '-', 'Color', colors(5,:),'LineWidth',3); 
plot(bestfit_amp.fitX, bestfit_amp.fitY_low, '-', 'Color', colors(3,:),'LineWidth',3); 
myerrorbar(orientationAll, highAveOcc, 'Color', colors(7,:), 'yError', highErrOcc,'Symbol','o','MarkerSize',10);
myerrorbar(orientationAll, medAveOcc, 'Color', colors(5,:), 'yError', medErrOcc,'Symbol','o','MarkerSize',10);
myerrorbar(orientationAll, lowAveOcc, 'Color', colors(3,:), 'yError', lowErrOcc ,'Symbol','o','MarkerSize',10);
set(gca,'FontSize', 13);
xlabel('Orientation Offset (Deg)');
ylabel('f1+f2  Amplitude (�V)');
title(sprintf('Amplitude Model Fobt = %0.4f', F_obt_fixed_joint));

subplot(2,2,3); hold on;
myerrorbar(orientationAll, noiseAveOcc2, 'Color', noiseColors, 'yError', noiseErrOcc2,'Symbol','o','MarkerSize',10);
plot(bestfit_full.fitX, bestfit_full.fitY_high2, '-', 'Color', colors(7,:),'LineWidth',3); 
plot(bestfit_full.fitX, bestfit_full.fitY_med2, '-', 'Color', colors(5,:),'LineWidth',3); 
plot(bestfit_full.fitX, bestfit_full.fitY_low2, '-', 'Color', colors(3,:),'LineWidth',3); 
myerrorbar(orientationAll, highAveOcc2, 'Color', colors(7,:), 'yError', highErrOcc2,'Symbol','o','MarkerSize',10);
myerrorbar(orientationAll, medAveOcc2, 'Color', colors(5,:), 'yError', medErrOcc2,'Symbol','o','MarkerSize',10);
myerrorbar(orientationAll, lowAveOcc2, 'Color', colors(3,:), 'yError', lowErrOcc2 ,'Symbol','o','MarkerSize',10);
set(gca,'FontSize', 13);
xlabel('Orientation Offset (Deg)');
ylabel('f2-f1  Amplitude (�V)');
title(sprintf('Full Model Fobt = %0.4f', F_obt_fixed_joint));

subplot(2,2,4); hold on;
myerrorbar(orientationAll, noiseAveOcc2, 'Color', noiseColors, 'yError', noiseErrOcc2,'Symbol','o','MarkerSize',10);
plot(bestfit_amp.fitX, bestfit_amp.fitY_high2, '-', 'Color', colors(7,:),'LineWidth',3); 
plot(bestfit_amp.fitX, bestfit_amp.fitY_med2, '-', 'Color', colors(5,:),'LineWidth',3); 
plot(bestfit_amp.fitX, bestfit_amp.fitY_low2, '-', 'Color', colors(3,:),'LineWidth',3); 
myerrorbar(orientationAll, highAveOcc2, 'Color', colors(7,:), 'yError', highErrOcc2,'Symbol','o','MarkerSize',10);
myerrorbar(orientationAll, medAveOcc2, 'Color', colors(5,:), 'yError', medErrOcc2,'Symbol','o','MarkerSize',10);
myerrorbar(orientationAll, lowAveOcc2, 'Color', colors(3,:), 'yError', lowErrOcc2 ,'Symbol','o','MarkerSize',10);
set(gca,'FontSize', 13);
xlabel('Orientation Offset (Deg)');
ylabel('f2-f1  Amplitude (�V)');
title(sprintf('Amplitude Model Fobt = %0.4f', F_obt_fixed_joint));
