for s = 1:nSubj
    thisFolder = dir(fullfile(dataDir, [subjectList{s},'*']));
    thisDir = fullfile(dataDir, thisFolder.name, matfileDir);
    subj_label = subjectList{s};
    
    [bestfit_full, bestfit_amp, F_joint_subj(s,1), AIC_joint_full_subj(s,1), AIC_joint_amp_subj(s,1)] = ...
        compareGaussianModels_IM(orientationAll, high_contrastsOcc(s,:), medium_contrastsOcc(s,:), low_contrastsOcc(s,:),...
        high_contrastsOcc2(s,:), medium_contrastsOcc2(s,:), low_contrastsOcc2(s,:));

end

figure;
F_sorted_joint = sort(F_joint_subj);
plot(F_sorted_joint,ones(length(F_sorted_joint),1)+(rand(length(F_sorted_joint),1)/10), 'o','color',[1 1 1], 'markerFaceColor',[0.2 0.2 0.2],'markerSize',15);
yaxis(0,2);
vline(2.361,'--');
xlabel('F-value');
set(gca,'ytick',[]);
title('joint');

AIC_diff_joint = AIC_joint_amp_subj - AIC_joint_full_subj;
figure;
AIC_sorted = sort(AIC_diff_joint);
plot(AIC_sorted, ones(length(AIC_sorted),1)+(rand(length(AIC_sorted),1)/10), 'o','color',[1 1 1], 'markerFaceColor',[0.2 0.2 0.2],'markerSize',15);
yaxis(0,2);
vline(0,'--');
xlabel('AIC diff (amp-full)');
set(gca,'ytick',[]);
title('joint');

figure; hold on;
plot(F_joint_subj, AIC_diff_joint,'o');
vline(2.361,'--');
hline(0,'--');
xlabel('F'); ylabel('AIC diff');



[bestfit_full, bestfit_amp, F_obt_joint,  AIC_joint_full, AIC_joint_amp] =compareGaussianModels_IM(orientationAll, highAveOcc, medAveOcc, lowAveOcc,...
    highAveOcc2, medAveOcc2, lowAveOcc2);

