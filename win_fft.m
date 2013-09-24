% FFT data for different time scale
% fft_scale: log2(windows)

function [fft_result] = win_fft(Data, fft_scale)

figure(fft_scale);

%Take out mean
Data = Data-mean(Data);

data_length = length(Data);
win_size = 2^fft_scale;
max_length = data_length - mod(data_length,win_size);

sample_num = max_length/win_size;

fft_result = zeros(1,win_size);

for i = 1 : sample_num
    fft_result = fft_result + fft(Data(1 + (i - 1) * win_size: i * win_size));
    % plot(abs(fftshift(fft_result)),'o');
end

plot(abs(fft_result),'x-');