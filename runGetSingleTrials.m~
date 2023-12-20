
dataDir = '~/data/orientationtuning';
for i = 1:nSubj
    thisFolder = dir(fullfile(dataDir, [subjectList{i},'*']));
    thisDir = fullfile(dataDir, thisFolder.name, matfileDir);
    [signalAll{i}, noiseLowerAll{i}, noiseUpperAll{i}, noiseMeanAll{i}] = getSingletrials(thisDir, subjectList{i}, taskList(i));
end

signal = cell(3,2);
noiseUpper = cell(3,2);
noiseLower = cell(3,2);
noiseMean = cell(3,2);
Contrasts = {'High','Med','Low'};
Orientations = [70, 90];
for i = 1:nSubj
    for iContrast = 1:3
        for iOrientation = 1:2
            signal{iContrast,iOrientation} = [signal{iContrast,iOrientation}; signalAll{i}{iContrast,iOrientation}];
            noiseUpper{iContrast,iOrientation} = [noiseUpper{iContrast,iOrientation}; noiseUpperAll{i}{iContrast,iOrientation}];
            noiseLower{iContrast,iOrientation} = [noiseLower{iContrast,iOrientation}; noiseLowerAll{i}{iContrast,iOrientation}];
            noiseMean{iContrast,iOrientation} = [noiseMean{iContrast,iOrientation}; noiseMeanAll{i}{iContrast,iOrientation}];
        end
    end
end
%%
x = 0:0.1:3;
y = 0:0.1:3;
figure;
for iContrast = 1:3
    for iOrientation = 1:2
        subplot(3,2,(iContrast-1)*2+iOrientation);
        plot(signal{iContrast,iOrientation}, noiseMean{iContrast,iOrientation}, 'o'); hold on;
        plot(x,y,'--','color',[0.5 0.5 0.5]);
        
        title(sprintf('%s %i',Contrasts{iContrast}, Orientations(iOrientation)));
        set(gca,'FontSize',15);
        xlabel('Signal Amplitude');
        ylabel('Noise Amplitude');
    end
end


%%
figure;
for iContrast = 1:3
    for iOrientation = 1:2
        subplot(3,2,(iContrast-1)*2+iOrientation);
        hSignal{iContrast,iOrientation} = histogram(signal{iContrast,iOrientation},'BinWidth',0.05); hold on;
        hLower{iContrast,iOrientation} = histogram(noiseLower{iContrast,iOrientation},'BinWidth',0.05);hold on;
        hUpper{iContrast,iOrientation} = histogram(noiseUpper{iContrast,iOrientation},'BinWidth',0.05);
        
        title(sprintf('%s %i',Contrasts{iContrast}, Orientations(iOrientation)));
        set(gca,'FontSize',15);
        xlabel('EEG Amplitude');
        ylabel('Count');
    end
end

for iContrast = 1:3
    for iOrientation = 1:2
        [h(iContrast,iOrientation), p(iContrast,iOrientation)] = kstest2(signal{iContrast,iOrientation}, noiseMean{iContrast,iOrientation});
    end
end