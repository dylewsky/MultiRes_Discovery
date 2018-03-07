clear variables; close all; clc

n = 1000;
fBounds = [1 100];
freqs = linspace(fBounds(1),fBounds(2),n);

amps = ones(1,n);
phases = 2*pi*rand(1,n);

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

save('uniform_data.mat','X','t','n','freqs');