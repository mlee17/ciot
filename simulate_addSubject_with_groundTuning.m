SNR_high0
secondSumAll = simulateTuning_tuningcurve;
tuningFunctions = [secondSumAll{1}; secondSumAll{2}; secondSumAll{3}];

nSub = length(SNR_high0);
nDraw = 100;
for iSelectSub = 1:nSub
    thisTuningFunctions = [];
    for iDraw = 1:nDraw
        thisDraw = randperm(nSub, iSelectSub);
        for iSub = 1:iSelectSub
            thisTuningFunctions(:,:,iDraw,iSub) = tuningFunctions * SNR_high0(thisDraw(iSub));
        end
        thisMeanTuningFunctions = mean(thisTuningFunctions,4);
    end
    meanTuningFunctions{iSelectSub} = mean(thisMeanTuningFunctions,3);
    meanTuningFunctions_se{iSelectSub} = std(thisMeanTuningFunctions,0,3) / sqrt(size(thisMeanTuningFunctions, 3));
end

orientationAll = [-90 -75 -45 -30 -15 -7 0 7 15 30 45 75 90];
colors = brewermap(7, 'BuPu');
noiseColors = brewermap(1,'Greens');
figure;
for iSelectSub = 1:10
    [bestfit_full{iSelectSub}, bestfit_amp{iSelectSub}, F_obt(iSelectSub), Fr(iSelectSub), AIC_full(iSelectSub), AIC_amp(iSelectSub)] = ...
        compareGaussianModels_fixedOffset(orientations, meanTuningFunctions{iSelectSub}(1,:), meanTuningFunctions{iSelectSub}(2,:), meanTuningFunctions{iSelectSub}(3,:));
    
    subplot(2,10,iSelectSub); hold on;
    plot(bestfit_amp{iSelectSub}.fitX, bestfit_amp{iSelectSub}.fitY_high, '-', 'Color', colors(7,:),'LineWidth',2); 
    plot(bestfit_amp{iSelectSub}.fitX, bestfit_amp{iSelectSub}.fitY_med, '-', 'Color', colors(5,:),'LineWidth',2); 
    plot(bestfit_amp{iSelectSub}.fitX, bestfit_amp{iSelectSub}.fitY_low, '-', 'Color', colors(3,:),'LineWidth',2); 
    myerrorbar(orientationAll, meanTuningFunctions{iSelectSub}(1,:), 'Color', colors(7,:), 'yError', meanTuningFunctions_se{iSelectSub}(1,:),'Symbol','o','MarkerSize',6);
    myerrorbar(orientationAll, meanTuningFunctions{iSelectSub}(2,:), 'Color', colors(5,:), 'yError', meanTuningFunctions_se{iSelectSub}(2,:),'Symbol','o','MarkerSize',6);
    myerrorbar(orientationAll, meanTuningFunctions{iSelectSub}(3,:), 'Color', colors(3,:), 'yError', meanTuningFunctions_se{iSelectSub}(3,:) ,'Symbol','o','MarkerSize',6);
    title(sprintf('Amp Model F=%0.4f', F_obt(iSelectSub)));
    box off;
    
    subplot(2,10,iSelectSub+10); hold on;
    plot(bestfit_full{iSelectSub}.fitX, bestfit_amp{iSelectSub}.fitY_high, '-', 'Color', colors(7,:),'LineWidth',2); 
    plot(bestfit_amp{iSelectSub}.fitX, bestfit_amp{iSelectSub}.fitY_med, '-', 'Color', colors(5,:),'LineWidth',2);
    plot(bestfit_amp{iSelectSub}.fitX, bestfit_amp{iSelectSub}.fitY_low, '-', 'Color', colors(3,:),'LineWidth',2); 
    myerrorbar(orientationAll, meanTuningFunctions{iSelectSub}(1,:), 'Color', colors(7,:), 'yError', meanTuningFunctions_se{iSelectSub}(1,:),'Symbol','o','MarkerSize',6);
    myerrorbar(orientationAll, meanTuningFunctions{iSelectSub}(2,:), 'Color', colors(5,:), 'yError', meanTuningFunctions_se{iSelectSub}(2,:),'Symbol','o','MarkerSize',6);
    myerrorbar(orientationAll, meanTuningFunctions{iSelectSub}(3,:), 'Color', colors(3,:), 'yError', meanTuningFunctions_se{iSelectSub}(3,:) ,'Symbol','o','MarkerSize',6);
   title(sprintf('Full Model F=%0.4f', F_obt(iSelectSub)));
   box off;
end

figure;
for iSelectSub = 11:20
    [bestfit_full{iSelectSub}, bestfit_amp{iSelectSub}, F_obt(iSelectSub), Fr(iSelectSub), AIC_full(iSelectSub), AIC_amp(iSelectSub)] = ...
        compareGaussianModels_fixedOffset(orientations, meanTuningFunctions{iSelectSub}(1,:), meanTuningFunctions{iSelectSub}(2,:), meanTuningFunctions{iSelectSub}(3,:));
    
    subplot(2,10,iSelectSub-10); hold on;
    plot(bestfit_amp{iSelectSub}.fitX, bestfit_amp{iSelectSub}.fitY_high, '-', 'Color', colors(7,:),'LineWidth',2); 
    plot(bestfit_amp{iSelectSub}.fitX, bestfit_amp{iSelectSub}.fitY_med, '-', 'Color', colors(5,:),'LineWidth',2); 
    plot(bestfit_amp{iSelectSub}.fitX, bestfit_amp{iSelectSub}.fitY_low, '-', 'Color', colors(3,:),'LineWidth',2); 
    myerrorbar(orientationAll, meanTuningFunctions{iSelectSub}(1,:), 'Color', colors(7,:), 'yError', meanTuningFunctions_se{iSelectSub}(1,:),'Symbol','o','MarkerSize',6);
    myerrorbar(orientationAll, meanTuningFunctions{iSelectSub}(2,:), 'Color', colors(5,:), 'yError', meanTuningFunctions_se{iSelectSub}(2,:),'Symbol','o','MarkerSize',6);
    myerrorbar(orientationAll, meanTuningFunctions{iSelectSub}(3,:), 'Color', colors(3,:), 'yError', meanTuningFunctions_se{iSelectSub}(3,:) ,'Symbol','o','MarkerSize',6);
    title(sprintf('Amp Model F=%0.4f', F_obt(iSelectSub)));
    box off;
    
    subplot(2,10,iSelectSub); hold on;
    plot(bestfit_full{iSelectSub}.fitX, bestfit_amp{iSelectSub}.fitY_high, '-', 'Color', colors(7,:),'LineWidth',2); 
    plot(bestfit_amp{iSelectSub}.fitX, bestfit_amp{iSelectSub}.fitY_med, '-', 'Color', colors(5,:),'LineWidth',2);
    plot(bestfit_amp{iSelectSub}.fitX, bestfit_amp{iSelectSub}.fitY_low, '-', 'Color', colors(3,:),'LineWidth',2); 
    myerrorbar(orientationAll, meanTuningFunctions{iSelectSub}(1,:), 'Color', colors(7,:), 'yError', meanTuningFunctions_se{iSelectSub}(1,:),'Symbol','o','MarkerSize',6);
    myerrorbar(orientationAll, meanTuningFunctions{iSelectSub}(2,:), 'Color', colors(5,:), 'yError', meanTuningFunctions_se{iSelectSub}(2,:),'Symbol','o','MarkerSize',6);
    myerrorbar(orientationAll, meanTuningFunctions{iSelectSub}(3,:), 'Color', colors(3,:), 'yError', meanTuningFunctions_se{iSelectSub}(3,:) ,'Symbol','o','MarkerSize',6);
   title(sprintf('Full Model F=%0.4f', F_obt(iSelectSub)));
   box off;
end