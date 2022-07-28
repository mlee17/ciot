function modelComparisonPlot(orientationAll, highAve, medAve, lowAve, noiseAve, highErr, medErr, lowErr, noiseErr, bestfit_full, bestfit_amp, F_obt, subjectID, freqLabel)

colors = brewermap(7, 'BuPu');
noiseColors = brewermap(1,'Greens');

figure; 
subplot(1,2,1); hold on;
myerrorbar(orientationAll, noiseAve, 'Color', noiseColors, 'yError', noiseErr,'Symbol','o','MarkerSize',10);
plot(bestfit_full.fitX, bestfit_full.fitY_high, '-', 'Color', colors(7,:),'LineWidth',3); 
plot(bestfit_full.fitX, bestfit_full.fitY_med, '-', 'Color', colors(5,:),'LineWidth',3); 
plot(bestfit_full.fitX, bestfit_full.fitY_low, '-', 'Color', colors(3,:),'LineWidth',3); 
myerrorbar(orientationAll, highAve, 'Color', colors(7,:), 'yError', highErr,'Symbol','o','MarkerSize',10);
myerrorbar(orientationAll, medAve, 'Color', colors(5,:), 'yError', medErr,'Symbol','o','MarkerSize',10);
myerrorbar(orientationAll, lowAve, 'Color', colors(3,:), 'yError', lowErr ,'Symbol','o','MarkerSize',10);
set(gca,'FontSize', 13);
xlabel('Orientation (Deg)');
ylabel([freqLabel,' EEG Power']);
title(sprintf('%s, Full Model Fobt = %0.4f', subjectID, F_obt));
subplot(1,2,2); hold on;
myerrorbar(orientationAll, noiseAve, 'Color', noiseColors, 'yError', noiseErr,'Symbol','o','MarkerSize',10);
plot(bestfit_amp.fitX, bestfit_amp.fitY_high, '-', 'Color', colors(7,:),'LineWidth',3); 
plot(bestfit_amp.fitX, bestfit_amp.fitY_med, '-', 'Color', colors(5,:),'LineWidth',3); 
plot(bestfit_amp.fitX, bestfit_amp.fitY_low, '-', 'Color', colors(3,:),'LineWidth',3); 
myerrorbar(orientationAll, highAve, 'Color', colors(7,:), 'yError', highErr,'Symbol','o','MarkerSize',10);
myerrorbar(orientationAll, medAve, 'Color', colors(5,:), 'yError', medErr,'Symbol','o','MarkerSize',10);
myerrorbar(orientationAll, lowAve, 'Color', colors(3,:), 'yError', lowErr ,'Symbol','o','MarkerSize',10);
set(gca,'FontSize', 13);
xlabel('Orientation (Deg)');
ylabel([freqLabel,' EEG Power']);
title(sprintf('%s, Amplitude Model Fobt = %0.4f', subjectID, F_obt));
