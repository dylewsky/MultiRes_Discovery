load('raw_data_5.mat')
load('raw_data_5_unmixed.mat');
load('mwDMD/mwDMD_sep_recon.mat');

% ground truth
gt_L = uv(:,1) * A(1,:);
gt_H = uv(:,2) * A(2,:);

% mwDMD
mw_L = xr_L;
mw_H = xr_H;
mw_t = tspan;

% Fourier
N = size(x,1);
dt = TimeSpan(2)-TimeSpan(1);
Fs = 1/dt; %sampling freq
freq = 0:Fs/N:Fs/2;

% % Compute power spectra 
x_power = zeros(length(freq),size(x,2));
for j = 1:size(x,2)
    xdft = fft(x(1:N,j));
    xdft = xdft(1:floor(N/2)+1);
    psdv = (1/(Fs*N)) * abs(xdft).^2;
    psdv(2:end-1) = 2*psdv(2:end-1);
    x_power(:,j) = psdv;
end
figure
plot(freq,x_power)
xlim([0 0.3])
xlabel('Frequency (Hz)')
ylabel('Power/Frequency (dB/Hz)')

freqThresh_L = 0.04;
freqThresh_H = 0.15;
steepness = 0.9;
% freq_bounds = [0.01 0.5];
% [ft_L, filt_L] = bandpass(x,[freq_bounds(1) freqThresh],Fs);
% [ft_H, filt_H] = bandpass(x,[freqThresh freq_bounds(2)],Fs);
[ft_L, filt_L] = lowpass(x,freqThresh_L,Fs,'Steepness',steepness);
[ft_H, filt_H] = highpass(x,freqThresh_H,Fs,'Steepness',steepness);
freqz(filt_L,1024,Fs)
freqz(filt_H,1024,Fs)

ft_L_power = zeros(length(freq),size(x,2));
ft_H_power = zeros(length(freq),size(x,2));
for j = 1:size(x,2)
    ft_L_dft = fft(ft_L(1:N,j));
    ft_L_dft = ft_L_dft(1:floor(N/2)+1);
    ft_L_psdv = (1/(Fs*N)) * abs(ft_L_dft).^2;
    ft_L_psdv(2:end-1) = 2*ft_L_psdv(2:end-1);
    ft_L_power(:,j) = ft_L_psdv;
    
    ft_H_dft = fft(ft_H(1:N,j));
    ft_H_dft = ft_H_dft(1:floor(N/2)+1);
    ft_H_psdv = (1/(Fs*N)) * abs(ft_H_dft).^2;
    ft_H_psdv(2:end-1) = 2*ft_H_psdv(2:end-1);
    ft_H_power(:,j) = ft_H_psdv;
end
figure
plot(freq,ft_L_power,'b')
hold on
plot(freq,ft_H_power,'r')
hold on
plot([freqThresh_L freqThresh_L],ylim,'b--')
hold on
plot([freqThresh_H freqThresh_H],ylim,'r--')
hold off
xlim([0 0.3])
xlabel('Frequency (Hz)')
ylabel('Power/Frequency (dB/Hz)')




tBounds = [mw_t(1) mw_t(end)];
lWidth = 1.5;

figure
subplot(2,2,1)
plot(TimeSpan,gt_L,'k','LineWidth',lWidth)
hold on
plot(mw_t,mw_L,'b','LineWidth',lWidth)
title('Slow-Scale: mwDMD vs Truth')
xlim(tBounds);

subplot(2,2,2)
plot(TimeSpan,gt_H,'k','LineWidth',lWidth)
hold on
plot(mw_t,mw_H,'b','LineWidth',lWidth)
title('Fast-Scale: mwDMD vs Truth')
xlim(tBounds);

subplot(2,2,3)
plot(TimeSpan,gt_L,'k','LineWidth',lWidth)
hold on
plot(TimeSpan,ft_L,'r','LineWidth',lWidth)
title('Slow-Scale: Fourier Filter vs Truth')
xlim(tBounds);

subplot(2,2,4)
plot(TimeSpan,gt_H,'k','LineWidth',lWidth)
hold on
plot(TimeSpan,ft_H,'r','LineWidth',lWidth)
title('Fast-Scale: Fourier Filter vs Truth')
xlim(tBounds);