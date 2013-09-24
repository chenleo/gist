% Binning data to I(t), correspond to the time_lag
% TimeHarp(s): 
% time_lag(ms): binning resolution
% route: which channel receive data
%   0:  for Donor
%   1:  for Acceptor
% Unit ms, first try time_lag = 0.1 ms

function [Donor, Acceptor, Total_DA] = binning(TimeHarp,route,time_lag)

Time_length = length(TimeHarp);         % length of the TimeHarp.
Time_max = int64(TimeHarp(Time_length));       % get the time length in s (experiment time)

% To_time_scale = nextpow2(Time_max*1000) - 1; 
To_use_ms = int64(Time_max*1000);             % Only use the time with the power of 2, unit, ms

Binning_length = To_use_ms / time_lag;  % Unit time_lag 0.1ms

% Initialize output
Donor = zeros(1,Binning_length);
Acceptor = zeros(1,Binning_length);
Total_DA = zeros(1,Binning_length);

%i = 1;
TimeHarp = TimeHarp * 1000/time_lag;
To_time_lag = To_use_ms/time_lag;
for i = 1:length(TimeHarp)
    DA_index = ceil(TimeHarp(i));
    if route(i) == 0
        Donor(DA_index) = Donor(DA_index) + 1;
    else
        Acceptor(DA_index) = Acceptor(DA_index) + 1;
    end
    Total_DA(DA_index) = Total_DA(DA_index) + 1;
    %i = i + 1;
    
end