% Test simulation for the contribution of spFRET from Triplet Event.

% Set Proximity FRET Ratio:
proximity = 0.6;

% Set Triplet effect:
% 0 for no triplet, 1 for all triplet;
D_tri = 1;
A_tri = 1;

% Set wave packet:
N = 1000;   % Number of sample time
W = 0.3;    % wave packet width, variance of wave packet
t = linspace(-1,1,N+1);
t = t(1:N); % Last point == First point

% Set wave packet of donor:
%% plot initial data;
figure(1);
a = subplot(2,2,1);
D = 1;
D_pack = D * exp(- t.^2/W);
D_trip = D_pack .* (1 - D_tri*rand(1,N));
plot(D_pack,'bo');
hold on;
plot(D_trip,'r-');
legend(a,'Prefect Donor wave packet','Donor wave packet with Triplet state');
title(a,'Simulate normalized wavepacket for Donor');
hold off;

% Set acceptor
% use definition of proximity = A/(D+A)
b = subplot(2,2,2);
if proximity == 1
    A = 0;
else 
    A = proximity * D/(1-proximity);
end
A_pack = A * exp(- t.^2/W);
A_trip = A_pack .* (1 - A_tri*rand(1,N));
plot(A_pack,'bo');
hold on;
plot(A_trip,'r-');
legend(b,'Prefect Acceptor wave packet','Acceptor wave packet with Triplet state');
title(b,'Simulate wavepacket for Acceptor, normalized related to Donor');
hold off;

% plot FRET value point(or hist) without triplet:
subplot(2,2,3);
% hist(A_pack./(D_pack+A_pack),1,'o');
plot(A_pack./(D_pack+A_pack),'o')

% plot FRET value point(or hist) with triplet:
subplot(2,2,4);
hist(A_trip./(D_trip+A_trip),100);
% plot(A_trip./(D_trip+A_trip),'o')

%% Plot proximity ratio after filtered 

figure(2);
filter_width = 2; % half width in frequency domain. If set to 1, only mean value will be left.
A_fft = fft(A_trip);
D_fft = fft(D_trip);

fft_filter = zeros(1,N);
fft_filter(1:filter_width) = 1;
fft_filter(N-filter_width+2:N) = 1;

% filter 
D_fft = D_fft.*fft_filter;
A_fft = A_fft.*fft_filter;

D_ifft = ifft(D_fft);
A_ifft = ifft(A_fft);
subplot(2,1,1);
result = A_ifft./(A_ifft + D_ifft);

plot(result(200:800),'o');
subplot(2,1,2);
hist(result(200:800),20,[0 1]);




