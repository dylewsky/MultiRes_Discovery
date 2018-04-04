clear variables; close all; clc

N = 100; %must be even
fBounds = [2 100];
fRes = 0.1;
freqs = fBounds(1):fRes:fBounds(2);

fPeaks = [24 61];
fSigma = [5 5];

amps = exp(-((freqs - fPeaks(1))/fSigma(1)).^2) + exp(-((freqs - fPeaks(2))/fSigma(2)).^2);
amps = amps/norm(amps);

tMax = 4 * 2*pi/min(freqs);
tStep = (1/32) * 2*pi/max(freqs);

t = 0:tStep:tMax;

X = zeros(N,length(t));

for n = 1:N
    phases = 2*pi*rand(size(amps));
    Y = repmat(amps.',1,length(t)).*sin(freqs.' * t + phases.');
    X(n,:) = sum(Y,1);
end

% 
% if 2*pi/min(freqs) > (t(end)-t(1))
%     disp('Not enough time for longest wavelength')
% end
% if 2*pi/max(freqs) < 16 * (t(2)-t(1))
%     disp('Fewer than 16 data points in one period of highest freq.')
% end

% X = repmat(amps.',1,length(t)).*sin(freqs.' * t + phases.');

save('bimodal_data.mat','X','t','N','freqs','fPeaks','fSigma');