% Calculate ACF by using fft for 1xN Vector data
% lag will be averaged

% Important: data should obey Wiener-Khinchin Theorem, (wide-sense-stationary random process)
% Fail to obey will lead to wrong result.

function acf = acf_fft(data)

data_mean = mean(data);

%calculate the fft for data
fft_data = fft(data);

%fft_data(1) is the 0 frequency, indicate the data mean.
acf = ifft(fft_data.*conj(fft_data)/fft_data(1))/data_mean -1 ;

%return the first half of acf
acf = acf(1:ceil(length(acf)/2));
