
dataDir = '~/data/orientationtuning';
for i = 1:nSubj
    thisFolder = dir(fullfile(dataDir, [subjectList{i},'*']));
    thisDir = fullfile(dataDir, thisFolder.name, matfileDir);
    [signalAll{i}, noiseLowerAll{i}, noiseUpperAll{i}] = getSingletrials(thisDir, subjectList{i}, taskList(i));
end

signal = cell(3,2);
noiseUpper = cell(3,2);
noiseLower = cell(3,2);
Contrasts = {'High','Med','Low'};
Orientations = [70, 90];
for i = 1:nSubj
    for iContrast = 1:3
        for iOrientation = 1:2
            signal{iContrast,iOrientation} = [signal{iContrast,iOrientation}; signalAll{i}{iContrast,iOrientation}];
            noiseUpper{iContrast,iOrientation} = [noiseUpper{iContrast,iOrientation}; noiseUpperAll{i}{iContrast,iOrientation}];
            noiseLower{iContrast,iOrientation} = [noiseLower{iContrast,iOrientation}; noiseLowerAll{i}{iContrast,iOrientation}];
        end
    end
end

figure;
for iContrast = 1:3
    for iOrientation = 1:2
        subplot(3,2,(iContrast-1)*2+iOrientation);
        hSignal{iContrast,iOrientation} = histogram(signal{iContrast,iOrientation},'BinWidth',0.05); hold on;
        hLower{iContrast,iOrientation} = histogram(noiseLower{iContrast,iOrientation},'BinWidth',0.05);hold on;
        hUpper{iContrast,iOrientation} = histogram(noiseUpper{iContrast,iOrientation},'BinWidth',0.05);
        
        title(sprintf('%s %i',Contrasts{iContrast}, Orientations(iOrientation)));
        set(gca,'FontSize',15);
    end
end