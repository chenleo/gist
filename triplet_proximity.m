% Test simulation for the contribution of spFRET from Triplet Event.

% Set Proximity FRET Ratio:
proximity = 0.5;

% Set Triplet effect:
% 0 for no triplet, 1 for all triplet;
D_tri = 0.7;
A_tri = 0.7;

% Set wave packet:
N = 1000;   % Number of sample time
W = 0.05;    % wave packet width, variance of wave packet
t = linspace(-1 , 1,N+1);
t = t(1:N); % Last point == First point

% Set wave packet of donor:
subplot(2,2,1);
D = 1;
D_pack = D * exp(- t.^2/W);
plot(D_pack);

% Set acceptor
% use definition of proximity = A/(D+A)
subplot(2,2,2);
if proximity == 1
    A = 0;
else 
    A = proximity * D/(1-proximity);
end
A_pack = A * exp(- t.^2/W);
plot(A_pack);

% plot FRET value point(or hist) without triplet:
subplot(2,2,3);
% hist(A_pack./(D_pack+A_pack),1,'o');
plot(A_pack./(D_pack+A_pack),'o')

% plot FRET value point(or hist) with triplet:
subplot(2,2,4);
A_trip = A_pack .* (1 - A_tri*rand(1,N));
D_trip = D_pack .* (1 - D_tri*rand(1,N));
hist(A_trip./(D_trip+A_trip));
% plot(A_trip./(D_trip+A_trip),'o')
