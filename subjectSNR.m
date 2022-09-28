subjectList = {'nl-1624','nl-jg','nl-1928','nl-0055','nl-1909','nl-1908','nl-0051','nl-0034','nl-myl','nl-rta','nl-0057','nl-2274','nl-2275',...
    'nl-2277','nl-2278','nl-2279','nl-2280','nl-2281','nl-2282', 'nl-2283','nl-2285', 'nl-2287', 'nl-2300','nl-1539','nl-2113','nl-2126','nl-2215',...
    'nl-2329','nl-2331','nl-2332','nl-2286'};%,'nl-2335','nl-2336','nl-2343','nl-2344','nl-2346','nl-2347','nl-2349','nl-2350','nl-2351','nl-2353', 'nl-2374',...
    %'nl-2378', 'nl-2379','nl-2381', 'nl-2383'};
subjectList = {'nl-2335','nl-2336','nl-2343','nl-2344','nl-2346','nl-2347','nl-2349','nl-2350','nl-2351','nl-2353', 'nl-2374',...
    'nl-2378', 'nl-2379','nl-2381', 'nl-2383'};
fileDir = '~/data/orientationRCA/snr';
dataDir = '~/data/orientationtuning';
matfileDir = 'Exp_MATL_HCN_128_Avg';
task = [zeros(31,1)];%;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1];
task = ones(15,1);
% for s = 1:length(subjectList)
%     loadDir = fullfile(fileDir, sprintf('%s_rcaSNR.mat',subjectList{s}));
%     load(loadDir);
%     rcaSNRAll(s,:) = rcaSNR;
% end
nSubj= length(subjectList);
meanSNR = [];
lowContrastSNR=[];
for i = 1:nSubj
    thisFolder = dir(fullfile(dataDir, [subjectList{i},'*']));
    thisDir = fullfile(dataDir, thisFolder.name, matfileDir);
    [meanSNR(i,:), allSNR_SumFreq(i,:),lowContrastSNR(i,:)] = computeSNR(thisDir, task(i));
end

% meanSNR(:,5:7) = rcaSNRAll;
% 
% meanSNR(:,8) = nanmean(meanSNR(:,1:7),2);
% meanSNR(:,9) = nanmean(meanSNR(:,1:4),2);
meanSNR(:,5) = nanmean(meanSNR(:,1:4),2);
meanSNR(length(subjectList)+1,:) = nanmean(meanSNR,1);

roundedSNR = round(meanSNR,2);
t = num2cell(roundedSNR);
t = cellfun(@num2str, t, 'UniformOutput', false);

x = repmat(1:size(meanSNR,2),nSubj+1,1);
y = repmat([1:nSubj+1]',1,size(meanSNR,2));

figure;
imagesc(meanSNR);
text(x(:),y(:), t, 'HorizontalAlignment','Center','Color',[0.8 0.8 0.8], 'FontSize',13);
set(gca, 'xTickLabel', {[],'2F1',[],'2F2',[],'1F1+1F2',[],'-1F+1F2',[],'Mean',[]});
set(gca,'yTick',1:nSubj+1);
set(gca, 'yTickLabel', [subjectList,'Mean']);
colorbar;
set(gca,'FontSize',13);
title('Occipital electrodes: High');
box off;
%%
% rank order
[orderedSNR, rankInd] = sort(meanSNR(1:nSubj,5), 'descend');
isOver2 = orderedSNR>2;

lowContrastSNR = lowContrastSNR(rankInd);
lowContrastSNR = lowContrastSNR(isOver2);
mean(lowContrastSNR)

allSNR_SumFreq_sorted = allSNR_SumFreq(rankInd, :);
allSNR_SumFreq_selected = allSNR_SumFreq_sorted(isOver2,:);
allSNR_SumFreq_mean = mean(allSNR_SumFreq_selected, 1);

highSNRSum = allSNR_SumFreq_mean(1:7);
medSNRSum = allSNR_SumFreq_mean(10:16);
lowSNRSum = allSNR_SumFreq_mean(19:25);

sumSNR = [highSNRSum; medSNRSum; lowSNRSum];

SNR_high0 = allSNR_SumFreq_selected(:,1);

%%
% correlation between raw IM and RCA IM
figure;
plot(meanSNR(1:10,3), meanSNR(1:10,6),'o', 'color',[1,1,1],'MarkerFaceColor',[0.4 0.4 0.8], 'MarkerSize',12);
title('(SNR 1F1+1F2) Raw vs RCA');
box off;
xlabel('Raw 1F1+1F2 SNR');
ylabel('RCA 1F1+1F2 SNR');
set(gca, 'FontSize',13);

figure;
plot(meanSNR(1:10,4), meanSNR(1:10,7),'o', 'color',[1,1,1],'MarkerFaceColor',[0.4 0.4 0.8], 'MarkerSize',12);
title('(SNR -1F1+1F2) Raw vs RCA');
box off;
xlabel('Raw -1F1+1F2 SNR');
ylabel('RCA -1F1+1F2 SNR');
set(gca, 'FontSize',13);

figure;
plot(meanSNR(1:10,5), meanSNR(1:10,8),'o', 'color',[1,1,1],'MarkerFaceColor',[0.4 0.4 0.8], 'MarkerSize',12);
title('(SNR -2F1+2F2) Raw vs RCA');
box off;
xlabel('Raw -2F1+2F2 SNR');
ylabel('RCA -2F1+2F2 SNR');
set(gca, 'FontSize',13);


figure;
plot(mean(meanSNR(1:10,3:5),2), mean(meanSNR(1:10,6:8),2),'o', 'color',[1,1,1],'MarkerFaceColor',[0.5 0.5 0.8], 'MarkerSize',12);
title('(SNR Mean IM) Raw vs RCA');
box off;
xlabel('Raw Mean IM SNR');
ylabel('RCA Mean IM SNR');
set(gca, 'FontSize',13);

% 
% %%%
% for i = 1:10
%     SNR76(i,:) = computeSNR;
% end
% 
% for i = 1:10
%     SNR76(i,7) = mean(SNR76(i,1:6));
% end
% 
% figure;
% imagesc(SNR76);
% set(gca, 'xTickLabel', {'2F1','2F2','F1+F2','-1F+1F2','-2F1+2F2','Mean'});
% colorbar;
% set(gca,'FontSize',13);
% title('Channel 76: High');

%%%
for i = 1:10
    meanSNRAll(i,:) = computeSNR2;
end

for i = 1:10
    meanSNRAll(i,7) = mean(meanSNRAll(i,1:6));
end

figure;
imagesc(meanSNRAll);
set(gca, 'xTickLabel', {'2F1','2F2','F1+F2','-1F+1F2','-2F1+2F2','Mean'});
colorbar;
set(gca,'FontSize',13);
title('All Contrasts');


%%%%
% corr between high and all contrasts
figure;
plot(meanSNR(:,7), meanSNRAll(:,7),'o', 'color',[1,1,1],'MarkerFaceColor',[0.4 0.4 0.8], 'MarkerSize',12);
title('High Contrast vs All Contrasts');
box off;
xlabel('High Contrast Mean SNR');
ylabel('All Contrasts Mean SNR');
set(gca, 'FontSize',13);

% corr between fundamental and Mean
figure;
plot(mean(meanSNRAll(:,1:2),2), meanSNRAll(:,7),'o', 'color',[1,1,1],'MarkerFaceColor',[0.5 0.5 0.8], 'MarkerSize',12);
title('Fundamental vs Mean SNR All Contrasts');
box off;
xlabel('Mean Fundamental All Contrasts');
ylabel('AMean SNR All Contrasts');
set(gca, 'FontSize',13);

% corr between fundamental and IM
figure;
plot(mean(meanSNRAll(:,1:2),2), mean(meanSNRAll(:,3:6),2),'o', 'color',[1,1,1],'MarkerFaceColor',[0.5 0.5 0.8], 'MarkerSize',12);
title('Fundamental vs IM All Contrasts');
box off;
xlabel('Mean Fundamental All Contrasts');
ylabel('Mean IM All Contrasts');
set(gca, 'FontSize',13);






