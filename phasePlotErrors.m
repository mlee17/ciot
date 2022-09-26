% function phasePlotErrors

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

nSubj= length(subjectList);

target_frequencies = [8,2];%[8, 2, 4]; %1f1+1f2 -1f1+1f2 -2f1+2f2
%orderedSNR, rankInd
orderedList = subjectList(rankInd);
for s = 1:length(orderedList)
    thisFolder = dir(fullfile(dataDir, [orderedList{s},'*']));
    thisDir = fullfile(dataDir, thisFolder.name, matfileDir);
        
    [output_cos{s}, output_sin{s}, output_amp{s}, x_est{s}, y_est{s}, V_indiv{s}, std_indiv{s}] = getPhases_trial(thisDir, target_frequencies);
    
end
% for s = 1:length(goodSubjList)
%     thisDir = fullfile(dataDir, goodSubjList{s}, matfileDir);
%     
%     [output_cos{s}, output_sin{s}, output_amp{s}, x_est{s}, y_est{s}, V_indiv{s}, std_indiv{s}] = getPhases_trial(thisDir, target_frequencies);
%     
% end
% 
% for s = 1:length(badSubjList)
%     thisDir = fullfile(dataDir, badSubjList{s}, matfileDir);
%     
%     [output_cos2{s}, output_sin2{s}, output_amp2{s}, x_est2{s}, y_est2{s}, V_indiv2{s}, std_indiv2{s}] = getPhases_trial(thisDir, target_frequencies);
%     
% end


goodColors = brewermap(13,'*BuGn');
goodColors = goodColors(1:9,:);
badColors = brewermap(4,'OrRd');

allColors = brewermap(27,'*BuGn');

% freqLabels = {'1F1+1F2', '-1F1+1F2', '-2F1+2F2'};
freqLabels = {'1F1+1F2', '-1F1+1F2'};
cycles = [1/8, 1/2];

figure;
for iFreq = 1:length(target_frequencies)
%     subplot(1,3,iFreq);
    subplot(1,2,iFreq);
    hold on;
    
    x_subj{iFreq} = []; y_subj{iFreq} = [];
    for s = 1:length(subjectList)
        plot([0 x_est{s}(iFreq)], [0 y_est{s}(iFreq)], 'color', allColors(s,:));
        hold on;
        circle(x_est{s}(iFreq), y_est{s}(iFreq), V_indiv{s}(iFreq), allColors(s,:));
        
        x_subj{iFreq}(s) = x_est{s}(iFreq);
        y_subj{iFreq}(s) = y_est{s}(iFreq);
    end
    
    meanX(iFreq) = mean(x_subj{iFreq});
    meanY(iFreq) = mean(y_subj{iFreq});
    plot([0 meanX(iFreq)], [0 meanY(iFreq)], 'color', [1, 0.2, 0.2], 'LineWidth',2);
    
    degree(iFreq) = atand(meanY(iFreq) / meanX(iFreq));
    radian(iFreq) = atan(meanY(iFreq) / meanX(iFreq));
    
    
%     for s = 1:length(goodSubjList)
%         plot([0 x_est{s}(iFreq)], [0 y_est{s}(iFreq)], 'color', goodColors(s,:));
%         hold on;
%         circle(x_est{s}(iFreq), y_est{s}(iFreq), V_indiv{s}(iFreq), goodColors(s,:));
%     end
%     for s = 1:length(badSubjList)
%         plot([0 x_est2{s}(iFreq)], [0 y_est2{s}(iFreq)], 'color', badColors(s,:));
%         hold on;
%         circle(x_est2{s}(iFreq), y_est2{s}(iFreq), V_indiv2{s}(iFreq), badColors(s,:));
%     end
    xlabel('cos');
    ylabel('sin');
%     yaxis(-1,1);
%     xaxis(-1,0.5);

    title(sprintf('%s, %0.4f',freqLabels{iFreq}, degree(iFreq)));
    set(gca,'FontSize',13);
    axis equal;
    hline(0);
    vline(0);
end

delay(1) = cycles(1) * ((180+degree(1))/360);
delay(2) = cycles(2) * (degree(2)/360);