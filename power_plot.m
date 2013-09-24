% Plot sum of energy from 0 point (mean), to max 

function [power_all] = power_plot(fft_result)
len_result = length(fft_result);
power_all = zeros(1,len_result/2);
power_map = fft_result.*conj(fft_result);   
power_sum = sum(power_map);
for i = 1:len_result/2
    power_all(i) = sum(power_map(1:i))/power_sum * 2;
end

figure(len_result);
plot(power_all,'x-');