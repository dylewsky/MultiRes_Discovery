clear variables; close all; clc

N = 100; %must be even
fBounds = [2 100];
freqs = linspace(fBounds(1),fBounds(2),N/2);
freqs = repelem(freqs,1,2);

% amps = ones(1,N);
amps = randn(1,N);

phases = 2*pi*rand(1,N);

tMax = 4 * 2*pi/min(freqs);
tStep = (1/32) * 2*pi/max(freqs);

t = 0:tStep:tMax;
% 
% if 2*pi/min(freqs) > (t(end)-t(1))
%     disp('Not enough time for longest wavelength')
% end
% if 2*pi/max(freqs) < 16 * (t(2)-t(1))
%     disp('Fewer than 16 data points in one period of highest freq.')
% end

X = repmat(amps.',1,length(t)).*sin(freqs.' * t + phases.');

save('pairwise_data.mat','X','t','N','freqs');