clear; close all; clc

%% Regular DMD

load('raw_data.mat');

[Phi, Omega, b] = run_dmd(x,TimeSpan);

xr = zeros(size(x));
for k = 1:size(x,1)
    xr = xr + Phi(:,k) * exp(Omega(k)*TimeSpan) * b(k);
%     figure
%     plot(TimeSpan,exp(Omega(k)*TimeSpan))
end

figure
scrollsubplot(2,1,1)
px = plot(TimeSpan,x);
title('Raw Data')
ylim(1.5*[min(min(x)) max(max(x))])
scrollsubplot(2,1,2)
pxr = plot(TimeSpan,real(xr));
title('DMD Reconstruction')
ylim(1.5*[min(min(x)) max(max(x))])

figure
plot(real(Omega),imag(Omega),'r*','LineWidth',2)
hold on
plot([-2 2],[0 0],'k--')
hold on
plot([0 0],[-2 2],'k--')
hold off
xlim([-2 2])
ylim([-2 2])
xlabel('Re(\lambda)')
ylabel('Im(\lambda)')

%% Multires DMD
% x = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16; -1 -2 -3 -4 -5 -6 -7 -8 -9 -10 -11 -12 -13 -14 -15 -16];
% x_PoT = x;

nLevels = 5;
nVars = size(x,1);
nSteps = 2^14;
x_PoT = x(:,1:nSteps); %truncate to a power of two
t_PoT = TimeSpan(:,1:nSteps);

mr_res = cell(nLevels,2^(nLevels-1));

nHold = 0;

for n = 1:nLevels
    nSplit = 2^(n-1);
    xL = reshape(x_PoT, nVars, nSteps/nSplit, nSplit);
    tL = reshape(t_PoT, 1, nSteps/nSplit, nSplit);
    
    for k = 1:nSplit
        xT = xL(:,:,k);
        tT = tL(:,:,k);
        mr_res{n,k}.x = xT;
        mr_res{n,k}.t = tT;
        [Phi, Omega, b] = run_dmd(xT,tT);
        mr_res{n,k}.Phi = Phi;
        mr_res{n,k}.Omega = Omega;
        mr_res{n,k}.b = b;
    end
%     nHold = input('Subtract off how many modes?')
end


%% Plot MultiRes Results
close all;
figure('units','pixels','Position',[0 0 1366 2*768])
plotDims = [3 4]; %rows, columns of plot grid on screen at a given time
for j = 1:nLevels
    nSplit = 2^(j-1);
    steps_per_window = nSteps/nSplit;
    om_spec = zeros(nVars,nSteps);
    b_spec = zeros(nVars,nSteps);

    scrollsubplot(plotDims(1),plotDims(2),[plotDims(2)*j-1, plotDims(2)*j]);
    plot(t_PoT,real(x_PoT),'k-') %plot ground truth
    xMax = max(max(abs(x_PoT)));
    
    ylim(1.5*[-xMax, xMax]);
    hold on
    
    for k = 1:nSplit
        Omega = mr_res{j,k}.Omega;
        Phi = mr_res{j,k}.Phi;
        b = mr_res{j,k}.b;
        t = mr_res{j,k}.t;
        tShift = t-t(1); %compute each segment of xr starting at "t = 0"
        
        om_sq = conj(Omega).*Omega;
        om_window_spec = repmat(om_sq, 1, steps_per_window);
        om_spec(:,(k-1)*steps_per_window+1:k*steps_per_window) = om_window_spec;
        
        scrollsubplot(plotDims(1),plotDims(2),plotDims(2)*j-3);
        plot(t,om_window_spec,'LineWidth',2);
%         title(['Frequency Spectrum for ' num2str(steps_per_window) '-Step Window']);
%         xlabel('t')
        xlim([t_PoT(1) t_PoT(end)]);
        ylabel('| \omega |^2')
        hold on

        if j ~=nLevels
            set(gca,'XTick',[])
        end
        
        b_sq = conj(b).*b;
        b_window_spec = repmat(b_sq, 1, steps_per_window);
        scrollsubplot(plotDims(1),plotDims(2),plotDims(2)*j-2);
        plot(t,b_window_spec,'LineWidth',2);
%         title(['Weights for ' num2str(steps_per_window) '-Step Window']);
%         xlabel('t')
        xlim([t_PoT(1) t_PoT(end)]);
        ylim auto;
        ylabel('|b|^2')
        hold on
        
        if j ~=nLevels
            set(gca,'XTick',[])
        end
        
        xr_window = zeros(nVars,steps_per_window);
        for q = 1:nVars
            xr_window = xr_window + Phi(:,q) * exp(Omega(q)*tShift) * b(q);
        end
        
        scrollsubplot(plotDims(1),plotDims(2),plotDims(2)*j-1:plotDims(2)*j);
        plot(t,real(xr_window),'LineWidth',2);
%         title(['Reconstruction for ' num2str(steps_per_window) '-Step Window']);
%         xlabel('t')
        xlim([t_PoT(1) t_PoT(end)]);
        ylabel('Re[x]')
        hold on

        if j ~=nLevels
            set(gca,'XTick',[])
        end
    end
    for k = 1:nSplit
        t = mr_res{j,k}.t;
        scrollsubplot(plotDims(1),plotDims(2),plotDims(2)*j-3);
    	plot([t(end) t(end)],get(gca, 'YLim'),'k:')
        hold on 
        scrollsubplot(plotDims(1),plotDims(2),plotDims(2)*j-2);
    	plot([t(end) t(end)],get(gca, 'YLim'),'k:')
        hold on 
        scrollsubplot(plotDims(1),plotDims(2),plotDims(2)*j-1:plotDims(2)*j);
    	plot([t(end) t(end)],get(gca, 'YLim'),'k:')
        hold on 
    end
end

%% Spectrogram
% for j = 1:nLevels
%     nSplit = 2^(j-1);
%     steps_per_window = nSteps/nSplit;
%     for k = 1:nSplit
%         
%     end
% end

sRate = 1/(t_PoT(2)-t_PoT(1));
figure
subplot(2,1,1)
plot(t_PoT,x_PoT,'LineWidth',2);
set(gca,'YTick',-4*pi:pi:2*pi);
grid on
% legend('x1','x2','x3','x4')
legend('a','b','r','\theta');
title('Time Series Data')
subplot(2,1,2)
for k = 1:nVars
%     scrollsubplot(2,2,k)
    
    k_PoT = fft(x_PoT(k,:));
    fSample = (0:size(k_PoT,2)-1)*sRate/size(k_PoT,2);
%     fShift = (-size(k_PoT,2)/2:size(k_PoT,2)/2-1)*(sRate/size(k_PoT,2));
%     kShift = fftshift(k_PoT);
    semilogx(fSample, abs(k_PoT),'LineWidth',2)
    xlabel('Frequency (Hz)')
    ylabel('Magnitude')
    xlim([fSample(1) fSample(end)]);
    
    hold on
%     title(['Fourier Spectrum of x' num2str(k)])

%     plot(fShift, abs(kShift))
end
title('Fourier Spectra')
% legend('x1','x2','x3','x4')
legend('a','b','r','\theta')
