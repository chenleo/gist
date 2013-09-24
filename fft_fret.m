% Script for Analyze FRET data
%% Data binning from TimeHarp and route
% seperate Donor, Acceptor, AD and binning 

time_lag = 0.4;              %ms
[Donor, Acceptor, Total_DA] = binning(TimeHarp,route,time_lag);




%% Calculate mean fft for different time window, and pick up the best one for filter
% for fft_scale = 4:6
%     [fft_result] = win_fft(Total_DA, fft_scale);
%     figure(fft_scale);
%     [power_all] = power_plot(fft_result);
% end

%% Calculate Proximity Ratio
filter_ratio = 0.32;
fft_scale = 6;
figure(1);
a = subplot(2,1,1);
hold on;
[Proximity, Burst_Begin, Burst_End, Burst_Size] = fft_filter(Total_DA, Total_DA, Acceptor, fft_scale, filter_ratio,10,3);
title(a,'Noremal FRET Hist');

b = subplot(2,1,2);
[Proximity2, Burst_Begin2, Burst_End2, Burst_Size2] = fft_filter(Acceptor, Total_DA, Acceptor, fft_scale, filter_ratio,10,3);

title(b,'Acceptor FRET Hist');
hold off;
median(Proximity2)
std(Proximity2)