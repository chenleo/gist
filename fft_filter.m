% FFT data for different time scale
% fft_scale: log2(windows) % for test, use 6: 2^6 = 64
% filter_ratio: sparcity of data to use % for 6: 0.2
function [Proximity, Burst_Begin, Burst_End, Burst_Size] = fft_filter(Ref,DA,Acceptor, fft_scale, filter_ratio, td,mw)

%figure(fft_scale);

% Take out mean
% Data = Data-mean(Data);

%% Initialize

% Total data length
data_length = length(DA);

% Window with
win_size = 2^fft_scale;

% Data is going to consider
max_length = data_length - mod(data_length,win_size);

% Max windows
sample_num = max_length/win_size;

% Set filter
fft_filter_point = floor(win_size/2*filter_ratio);
fft_filter = zeros(1,win_size);

% Initial filter
fft_filter(1:fft_filter_point) = 1; 
fft_filter(win_size + 2 - fft_filter_point : win_size) = 1;

% Initail Output
Proximity = 0;
Burst_Begin = 0;
Burst_End = 0;

% Count burst
Burst_Count = 0;

%% main loop to select burst
% Burst selection procedure will be taken in filtered data
% While Proximity ratio will be calculate using the raw_data

threshold = td;
max_width = mw;
for i = 1 : sample_num
    loop_start = 1 + (i - 1) * win_size;
    loop_end = i* win_size;
    fft_result = fft(Ref(loop_start: loop_end));
    fft_filtered = fft_result .* fft_filter;
    
    % data to evaluate
    data_filtered = ifft(fft_filtered);
    [max_value, max_ind] = max(data_filtered);
    
    if max_value > threshold && max_ind - fft_scale/2 > 0 && max_ind + fft_scale/2 <= win_size
        Burst_Count = Burst_Count + 1;
        Burst_Begin(Burst_Count) = max_ind - max_width + (i - 1) * win_size;
        Burst_End(Burst_Count) = max_ind + max_width + (i - 1) * win_size;
        Burst_Size(Burst_Count) = sum(DA(Burst_Begin(Burst_Count):Burst_End(Burst_Count)));
        
        % P = A/(A+D)
        Proximity(Burst_Count) = sum(Acceptor(Burst_Begin(Burst_Count):Burst_End(Burst_Count)))/sum(DA(Burst_Begin(Burst_Count):Burst_End(Burst_Count)));
        
    end
    
%command out for debugging 
%     plot(Ref(1 + (i - 1) * win_size: i * win_size),'o-');
%     hold on
%     plot(DA(1 + (i - 1) * win_size: i * win_size),'ko-');
%     plot(data_filtered, 'rx-');
%     plot(ones(1,win_size) * threshold, 'g.');
%     % plot(abs(fftshift(fft_result)),'o');
%     hold off;
end
hist(Proximity,50);