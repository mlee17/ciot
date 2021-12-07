occipital_electrodes = [69,70,73,74,71,75,81,76,82,83,84,89,88];

for s = 1:nSubj
    thisFolder = dir(fullfile(dataDir, [subjectList{s},'*']));
    thisDir = fullfile(dataDir, thisFolder.name, matfileDir);
    subj_label = subjectList{s};
    
    output = loadAxx(2, thisDir);
    
    thisSpec = output.freq_ampl(:, occipital_electrodes, 1); % only condition 1
    thisSpec = mean(thisSpec, 2);
    
    specAll(:,s) = thisSpec;

    thisCos = output.cos(:, occipital_electrodes, 1); % only condition 1
    thisCos = mean(thisCos, 2);
    thisSin = output.sin(:, occipital_electrodes, 1);
    thisSin = mean(thisSin, 2);
    
    cosAll(:,s) = thisCos;
    sinAll(:,s) = thisSin;
%     
end

specMean = mean(specAll, 2);
frequencies = 0:0.5:50;

cosMean = mean(cosAll, 2);
sinMean = mean(sinAll, 2);
specMean = sqrt(cosMean.^2 + sinMean.^2);
frequencies = 0:0.5:50;
%%%
thisDir = '~/proj/2021_orientationSingleFrequency/Project_20211201_1747/Exp_MATL_HCN_128_Avg';
output = loadAxx(2, thisDir);
specMean = output.freq_ampl(:, occipital_electrodes, 2);
%%%
figure;
stem(frequencies, specMean(1:101),'Color',[0.5 0.5 0.5], 'Marker','none', 'LineWidth', 1.5);
hold on;
f1 = 3:3:48;
f1_index = f1/0.5 + 1;
stem(f1, specMean(f1_index), 'Color', [0.2 0.7 0.3], 'Marker', 'none', 'LineWidth', 3);

f2 = 5:5:50;
f2_index = f2/0.5 + 1;
stem(f2, specMean(f2_index), 'Color', [0.2 0.3 0.7], 'Marker', 'none', 'LineWidth', 3);

IM = [1,2,4,8,16];
IM_index = IM/0.5 + 1;
stem(IM, specMean(IM_index), 'Color', [0.9 0.4 0.8], 'Marker', 'none', 'LineWidth', 3);

box off;
legend({'','F1','F2','IM'});
xlabel('Frequency (Hz)');
ylabel('Amplitude');


