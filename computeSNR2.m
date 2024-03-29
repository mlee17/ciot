function meanSNR = computeSNR2%(cond, conditions_to_visualize, groups_to_visualize, target_freq, file_name)
% All contrasts
cond = 2;
conditions_to_visualize = '1-27';
groups_to_visualize = {'8-9','17-18','26-27'};%{'1-7', '10-16', '19-25'};
target_freq1= 6; % 2F1
target_freq2= 10; % 2F2
occipital_electrodes = [69,70,73,74,71,75,81,76,82,83,84,89,88];
freq = [6 10 8 2 4]; % 2f1 2f2 f1+f2 -f1+f2 -2f1+2f2
% ex) significantSNRPlot (1, '1-27', {'1-7', '12-15', '19-26'}, 10)
% this function plots the orientation versus frequency amplitudes of EEG data 
% for electrodes with a significant signal to noise ratio
% split into different specified contrasting levels obtained from PowerDiva up to 50 Hz of frequency.

% ***********
% Variables
% ***********

% cond relates to the sort of data you would like to analyze. Your options are: 
%  1) Axx_trial 
%  2) Axx 
%  3) Raw EEG

% conditions to visualize is the conditions within the directory you would
% like to visualize. This is variable, depending on your data. 

% groups to visualize are the conditions that you want to group together for contrasting levels
% ex: ['1-7', '10-16', '19-25'] corresponds to low, medium and high contrasting

% target_freq represents the frequency that you want to graph the amplitudes for

% file_name represents the name of the file to save the plot to. If not provided, the default 
% file name of significantSNRPlot will be used

%***********************
% running this function
% **********************

% significantSNRPlot (3, '1-27', {'1-7', '10-16', '19-25'}, 10, 'mySNRPlot')

% the above command filters (SNR ratio > 4), averages and then visualizes the electrode data for
% Raw EEG data (cond = 3) for
% conditions 1 to 27 (conditions_to_visualize)
% where the frequency is equal to 10 (target_freq)
% with the low contrast set from '1-7', medium contrast '10-16', and high contrast '19-25'

% The plot will be saved in a file called 'mySNRPlot'

% --------------------------------------------------------------------

% this block checks the conditions_to_visualize entry and splits it into a
% matrix of condition numbers
if exist('conditions_to_visualize', 'var') == 0 
    disp('please enter a valid condition 1-9 ');
    return
elseif ~isempty(strfind(conditions_to_visualize, ','))
    
    conditions_to_visualize = str2double(strsplit(conditions_to_visualize,','));
    
elseif ~isempty(strfind(conditions_to_visualize, '-'))
    
    conditions_to_visualize = str2double(strsplit(conditions_to_visualize,'-'));
    conditions_to_visualize = conditions_to_visualize(1):conditions_to_visualize(2);
else 
    conditions_to_visualize = str2double(conditions_to_visualize);
end

% this block checks if the groups_to_visualize variable exists
% and then converts the string input into the corresponding range and groups them into parsed_group
% ex: if groups_to_visualize is {'1-7', '12-15', '19-26'}
% then parsed_group = [1, 2, 3, 4, 5, 6, 7; 12, 13, 14, 15; 19, 20, 21, 22, 23, 24, 25, 26]
parsed_group = [];

if exist('groups_to_visualize', 'var') == 0
    disp('please enter an array with groups to visualize')
    return
else 
    parsed_group = [];
    jj = 1;
    for groups = groups_to_visualize
        % each element in groups_to_visualize must contain a dash (ex: '1-9')
        groups = str2double(strsplit(groups{1},'-'));
        groups = groups(1):groups(2);
        parsed_group = [parsed_group; groups];
        jj = jj + 1;
    end
end

% Check for file_name or assign default file name
if exist('file_name', 'var') == 0
    file_name = 'significantSNRPlot';
end

% Converts the condition to a number
if lower(cond) == 'raw'
    cond = 3;
elseif lower(cond) == 'trial'
    cond = 1;
elseif lower(cond) == 'axx'
    cond = 2;
end

% Reads Power Diva Data
if cond == 3
    data = readPowerDiva(cond);
else
    [~, ampl] = readPowerDiva(cond);
end

% Verifies the correctness of data with the conditions
if exist('data','var') == 1
    o = size(data);
    if o(end) < max(conditions_to_visualize) % checks to make sure the number of conditions the user enters matchthe data in the directory
        disp('At least one of your selected conditions exceeds the number of conditions in this directory.')
        disp('Please revise the number of conditions.');
        disp('Quitting the program!');
        return
    end
else
    o = size(ampl);
    if o(end) < max(conditions_to_visualize) % checks to make sure the number of conditions the user enters matchthe data in the directory
    disp('At least one of your selected conditions exceeds the number of conditions in this directory.')
    disp('Please revise the number of conditions.');
    disp('Quitting the program!');
    return
    end
end


% this if block determines the size of the subplot plot
if length(conditions_to_visualize) == 1
    [subplot_x, subplot_y] = subplot_num_gen(1);
    conditions_to_visualize = [conditions_to_visualize;conditions_to_visualize]; %this ensures that the loop below will work
else
    [subplot_x, subplot_y] = subplot_num_gen(length(conditions_to_visualize)); % this finds the
    % optimal dimensions for the subplot
end


%*************************************************************
%*************************************************************

% this section is for reading the raw EEG data

%*************************************************************
%*************************************************************

% Extract the ranges for the different contrast values
low_contrast = parsed_group(3,:);

medium_contrast = parsed_group(2,:);

high_contrast = parsed_group(1,:);

orientation = [-90,-70,-45,-30,-15,-7, 0, 7, 15, 30, 45, 70, 90]; % negative values included for reflection
signal_index = 1;
all_snrs = 0;

for channel_to_visualize = 1 : 128

    if cond == 3
        a = squeeze(data(:,2:size(data,2) - 1, channel_to_visualize, :,:)); % excluding the first and the last epochs from the analysis
        % a = squeeze(data(:,:, channel_to_visualize, :,:));
        new_dimensions = [size(data,1)*2, (size(data,2)-2)/2, size(data,4), size(data,5)]; % reshaping it because we want to
        % use 2 second epochs for plotting frequency ampl. This increases the
        % resolution to 0.5 hz as opposed to 1 hz.
        a = reshape(a, new_dimensions);
        a = squeeze(mean(mean(a,2),3)); % averaging epochs
        jj = 1; % loop counter
        
        if exist('ampls', 'var') ==1
            clear ampls
        end
        
    %     ampls = zeros(str2double(conditions_to_visualize(2))...
    %         - str2double(conditions_to_visualize(1))+1,2);
        for i = 1: length(conditions_to_visualize)
            ampls(i,:) = freq_ampl(a(:,i), 420, 100);
        end
        
        ampls = max(max(ampls)); % this variable is used to normalize the scale of
        % y-axis across all graphs
        
        for c = conditions_to_visualize % looping over conditions 
            x = 0:1000/size(data,1):1000;
            x = x(1:size(data,1)); %timestamp for graphing
            % Get all the amplitudes corresponding to target_freq within the jj-th condition
            % and sampling frequency of 420
            amplitude = get_frequency(a(:, conditions_to_visualize(jj)), 420, target_freq);
            average_amplitude = get_avg_frequency(a(:, conditions_to_visualize(jj)), 420, target_freq);
            % condition_amplitude stores the amplitudes to be graphed (from 1 to number of conditions)
            condition_amplitude(jj) = amplitude;
            average_condition_amplitude(jj) = average_amplitude;
            jj = jj+1;
        end


    %*************************************************************
    %*************************************************************

    % this section is for plotting freq ampl of processed Axx files w.o trial data
    % The looping process is similar to condition 3, where we get the corresponding
    % amplitudes for a target_freq, but with the data coming from Axx files instead.

    %*************************************************************
    %*************************************************************
        
    elseif cond == 2
        jj = 1;
        a = squeeze(ampl (:, channel_to_visualize, :));
        ampl_scale = a(:, conditions_to_visualize);
        ampl_scale = max(max(ampl_scale));

        for c = conditions_to_visualize
            for iFreq = 1:length(freq)
                thisAmp = axx_get_frequency(a(:,conditions_to_visualize(jj)), freq(iFreq));
                allAmp(iFreq,jj) = thisAmp;
                [left right] = axx_get_sidebands(a(:, conditions_to_visualize(jj)), freq(iFreq));
                thisSidebands = (left+right)/2;
                allSidebands(iFreq, jj) = thisSidebands;
            end
                jj = jj + 1;
%                 
%             amplitude = axx_get_frequency(a(:,conditions_to_visualize(jj)), target_freq1);
%             condition_amplitude1(jj) = amplitude;
%             amplitude = axx_get_frequency(a(:,conditions_to_visualize(jj)), target_freq2);
%             condition_amplitude2(jj) = amplitude;
%             conditions(jj) = jj;
%             [left right] = axx_get_sidebands(a(:, conditions_to_visualize(jj)), target_freq1);
% %             average_condition_amplitude(jj) = average_amplitude;
%             sidebands1(jj,:) = [left right];
%             [left right] = axx_get_sidebands(a(:, conditions_to_visualize(jj)), target_freq2);
%             sidebands2(jj,:) = [left right];
%             
%             jj = jj + 1;
        end

    %*************************************************************
    %*************************************************************

    % this section is for reading the processed Axx files WITH trial data
    % The looping process is similar to condition 3, where we get the corresponding
    % amplitudes for a target_freq, but with the data coming from Axx files with trial data instead.

    %*************************************************************
    %*************************************************************

        
    else
        jj = 1;
        a = squeeze(mean(ampl,3));
        a = squeeze(a (:, channel_to_visualize, :));
        
        ampl_scale = a(:, conditions_to_visualize);
        ampl_scale = max(max(ampl_scale));
        
        for c = conditions_to_visualize % looping over conditions
            
            amplitude = axx_get_frequency(a(:,conditions_to_visualize(jj)), target_freq);
            condition_amplitude(jj) = amplitude;
            average_amplitude = axx_get_avg_frequency(a(:, conditions_to_visualize(jj)), target_freq);
            average_condition_amplitude(jj) = average_amplitude;
            conditions(jj) = jj;
            jj = jj +1;
        end



    end

    % condition_amplitude is the corresponding amplitude of each contrast variable
    % We reflect the corresponding amplitudes to account for the negative degrees (-90 to 0)

    % Compute the average of low, medium and high for the average_condition_amplitudes
%     average_condition_amplitude = split(average_condition_amplitude, low_contrast, medium_contrast, high_contrast);
% 
%     % Compute max value among three contrasts and compare it to average
%     low_contrasts = reflect(condition_amplitude(low_contrast));
%     medium_contrasts = reflect(condition_amplitude(medium_contrast));
%     high_contrasts = reflect(condition_amplitude(high_contrast));
%     average_contrasts = reflect(average_condition_amplitude);
% 
%     % Compute the signal to noise ratio
%     snr_ratio = find_SNR(orientation, low_contrasts, medium_contrasts, high_contrasts, average_contrasts);
%   

    %2F1 -> high contrast only Cond 8
    allAmp(1,8);
    allSidebands(1,8);
    %2F2 -> high contrast only Cond 9
    allAmp(2,9);
    allSidebands(2,9);
    %-1F1+1F2 -> high contrast only Cond 1
    allAmp(3,1);
    allSidebands(3,1);
    %-2F1+2F2 -> high contrast only Cond 1
    allAmp(4,1);
    allSidebands(4,1);
    %F1+F2 -> high contrast only Cond 1
    allAmp(5,1);
    allSidebands(5,1);
    %2F1+2F2 -> high contrast only Cond 1
    allAmp(6,1);
    allSidebands(6,1);
    
    %high
    thisChannelAmp(1,:) = [allAmp(1,8), allAmp(2,9), allAmp(3,1), allAmp(4,1), allAmp(5,1), allAmp(6,1)];
    %med
    thisChannelAmp(2,:) = [allAmp(1,17), allAmp(2,18), allAmp(3,10), allAmp(4,10), allAmp(5,10), allAmp(6,10)];
    %low
    thisChannelAmp(3,:) = [allAmp(1,26), allAmp(2,27), allAmp(3,19), allAmp(4,19), allAmp(5,19), allAmp(6,19)];
    
    thisChannelSide(1,:) = [allSidebands(1,8), allSidebands(2,9), allSidebands(3,1), allSidebands(4,1), allSidebands(5,1), allSidebands(6,1)];
    thisChannelSide(2,:) = [allSidebands(1,17), allSidebands(2,18), allSidebands(3,10), allSidebands(4,10), allSidebands(5,10), allSidebands(6,10)];
    thisChannelSide(3,:) = [allSidebands(1,26), allSidebands(2,17), allSidebands(3,19), allSidebands(4,19), allSidebands(5,19), allSidebands(6,19)];

    thisChannelAmp = mean(thisChannelAmp, 1);
    thisChannelSide = mean(thisChannelSide,1);
%     high1 = [sidebands1(high_contrast(1), 1), condition_amplitude1(high_contrast(1)), sidebands2(high_contrast(1),2)];
%     high2 = [sidebands2(high_contrast(2), 1), condition_amplitude1(high_contrast(2)), sidebands2(high_contrast(2),2)];
%     
%     med1 = [sidebands1(medium_contrast(1), 1), condition_amplitude1(medium_contrast(1)), sidebands2(medium_contrast(1),2)];
%     med2 = [sidebands2(medium_contrast(2), 1), condition_amplitude1(medium_contrast(2)), sidebands2(medium_contrast(2),2)];
%     
%     low1 = [sidebands1(low_contrast(1), 1), condition_amplitude1(low_contrast(1)), sidebands2(low_contrast(1),2)];
%     low2 = [sidebands2(low_contrast(2), 1), condition_amplitude1(low_contrast(2)), sidebands2(low_contrast(2),2)];
%     

%     if snr_ratio > 4
    if any(channel_to_visualize == occipital_electrodes)
        % Include this electrode's data if the SNR ratio is high enough
%         significant_low_contrasts(signal_index, :) = low_contrasts;
%         significant_medium_contrasts(signal_index, :) = medium_contrasts;
%         significant_high_contrasts(signal_index, :) = high_contrasts;
%         significant_average_contrasts(signal_index, :) = average_contrasts;
        occ_Amp(signal_index,:) = thisChannelAmp;
        occ_Side(signal_index,:) = thisChannelSide;
        
%         occ_low1(signal_index,:) = low1;
%         occ_low2(signal_index,:) = low2;
%         occ_med1(signal_index,:) = med1;
%         occ_med2(signal_index,:) = med2;
%         occ_high1(signal_index,:) = high1;
%         occ_high2(signal_index,:) = high2;

        signal_index = signal_index + 1;
%         all_snrs = all_snrs + snr_ratio;
    end

end % end of for loop

allSNR = occ_Amp./occ_Side;
meanSNR = mean(allSNR, 1);
% 
% colors = brewermap(7, 'BuPu');
% noiseColors = brewermap(1,'Greens');
% 
% 
% figure;
% subplot(3,1,1);
% bar([mean(occ_high1);mean(occ_high2)], 'FaceColor', colors(7,:));
% set(gca,'xTickLabel',{'2F1','2F2'})
% ylabel('high')
% % bar(mean(occ_high1));
% % subplot(3,2,2);
% % bar(mean(occ_high2));
% 
% subplot(3,1,2);
% bar([mean(occ_med1);mean(occ_med2)], 'FaceColor', colors(5,:));
% set(gca,'xTickLabel',{'2F1','2F2'})
% ylabel('med')
% % bar(mean(occ_med1));
% % subplot(3,2,4);
% % bar(mean(occ_med2));
% 
% subplot(3,1,3);
% bar([mean(occ_low1);mean(occ_low2)], 'FaceColor', colors(3,:));
% set(gca,'xTickLabel',{'2F1','2F2'})
% ylabel('low')
% % bar(mean(occ_low1));
% % subplot(3,2,6);
% % bar(mean(occ_low2));


% 
% if(signal_index == 1)
%     disp('No electrode with significant signal-to-noise ratio')
% else
%     % Compute the average of all electrode data
%     mean_low = mean(significant_low_contrasts);
%     mean_medium = mean(significant_medium_contrasts);
%     mean_high = mean(significant_high_contrasts);
%     mean_average = mean(significant_average_contrasts);
%     miny = min([mean_low, mean_medium, mean_high, mean_average]);
%     maxy = max([mean_low, mean_medium, mean_high, mean_average]);
% 
%     % snr_ratio of graph
%     snr_ratio = find_SNR(orientation, mean_low, mean_medium, mean_high, mean_average);
%     % average of all snrs
%     all_snr_ratio = all_snrs / (signal_index + 1);
% 
% colors = brewermap(7, 'BuPu');
% noiseColors = brewermap(1,'Greens');
% 
%     fig = figure();
%     
%     hold on
%     
%     plot(orientation, mean_high,'Color', colors(7,:))
%     plot(orientation, mean_medium,'Color', colors(5,:))
%     plot(orientation, mean_low,'Color', colors(3,:))
%     plot(orientation, mean_average,'Color', noiseColors)
%     ylim([miny maxy])
%     title('EEG Power vs Orientation')
%     xlabel('Orientation (Degrees)')
%     ylabel('F1+F2 EEG power (�V)')
%     legend('High Contrast','Medium Contrast','Low Contrast', 'NoiseAverage')
%     text_box_dim = [.2 .5 .3 .3];
%     str = sprintf('Ratio of max heights is %d\n Calculated ratio is %d', snr_ratio, all_snr_ratio);
%     annotation('textbox',text_box_dim,'String',str,'FitBoxToText','on');
% %     print(fig,file_name,'-dpng', '-r300')
% end
% 
% 
end

% This function 'reflects' an array. 
% Used to prepare the amplitudes to be graphed with negative orientations.
% ex: reflect([1, 2, 3]) => [3, 2, 1, 2, 3]
function y = reflect(a)
    y = [fliplr(a(2:end)), a];
end

% Splits an array of frequencies into three contrasts and averages them
% a is the array, length if the length of each contrast
% ex: 
% low_contrast = [1,2]
% medium_contrast = [3,4]
% high_contrast = [5,6]
% split([1, 2, 5, 6, 3, 10], low_contrast, medium_contrast, high_contrast) => [3, 6]
% function assumes that there are three groupings
function x = split(a, low_contrast, medium_contrast, high_contrast)
    low = a(low_contrast);
    med = a(medium_contrast);
    high = a(high_contrast);
    x = (low + med + high) / 3;
end