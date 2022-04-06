figure;
for s = 1:nSubj
    subplot(5,4,s); hold on;
    plot(orientationAll, average_contrastsOcc(s,:), '-', 'Color', noiseColors);
    plot(orientationAll, high_contrastsOcc(s,:),'-', 'Color', colors(7,:));
    plot(orientationAll, medium_contrastsOcc(s,:),'-', 'Color', colors(5,:));
    plot(orientationAll, low_contrastsOcc(s,:),'-', 'Color', colors(3,:));
    
    myerrorbar(orientationAll, average_contrastsOcc(s,:), 'Color', noiseColors, 'Symbol','o','MarkerSize',6);
    myerrorbar(orientationAll, high_contrastsOcc(s,:), 'Color', colors(7,:), 'Symbol','o','MarkerSize',6);
    myerrorbar(orientationAll, medium_contrastsOcc(s,:), 'Color', colors(5,:), 'Symbol','o','MarkerSize',6);
    myerrorbar(orientationAll, low_contrastsOcc(s,:), 'Color', colors(3,:), 'Symbol','o','MarkerSize',6);

        subj_label = subjectList{s};
    title(subj_label);
%     xlabel('Orientation Offset (Deg)');
%     ylabel('f2-f1 Amplitude (µV)');
end
figure;
for s = 1:nSubj
    subplot(5,4,s); hold on;
    plot(orientationAll, average_contrastsOcc2(s,:), '-', 'Color', noiseColors);
    plot(orientationAll, high_contrastsOcc2(s,:),'-', 'Color', colors(7,:));
    plot(orientationAll, medium_contrastsOcc2(s,:),'-', 'Color', colors(5,:));
    plot(orientationAll, low_contrastsOcc2(s,:),'-', 'Color', colors(3,:));
    
    myerrorbar(orientationAll, average_contrastsOcc2(s,:), 'Color', noiseColors, 'Symbol','o','MarkerSize',6);
    myerrorbar(orientationAll, high_contrastsOcc2(s,:), 'Color', colors(7,:), 'Symbol','o','MarkerSize',6);
    myerrorbar(orientationAll, medium_contrastsOcc2(s,:), 'Color', colors(5,:), 'Symbol','o','MarkerSize',6);
    myerrorbar(orientationAll, low_contrastsOcc2(s,:), 'Color', colors(3,:), 'Symbol','o','MarkerSize',6);

    subj_label = subjectList{s};
    title(subj_label);
%     xlabel('Orientation Offset (Deg)');
%     ylabel('f2-f1 Amplitude (µV)');
end


