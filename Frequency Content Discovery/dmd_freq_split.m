clear variables; close all; clc

%% Initialize
N = 100; %must be even
fBounds = [2 100];
freqs_unsplit = linspace(fBounds(1),fBounds(2),N/2);
freqs_unsplit = repelem(freqs_unsplit,1,2);

% amps = ones(1,N);
amps = rand(1,N);

phases = 2*pi*rand(1,N);

tMax = 4 * 2*pi/min(freqs_unsplit);
tStep = (1/32) * 2*pi/max(freqs_unsplit);

t = 0:tStep:tMax;

%define which DMD function to use
dmd = @(x,t,r) reg_dmd(x,t,r);

splits = 0:0.1:0.5;
wSize = 100;
stepSize = 10;
nSteps = floor((length(t) - wSize)/stepSize);

vs_res = cell(length(splits),nSteps);
vs_freqs = cell(length(splits),1);
    
r = N;
%% Generate Split-Freq Data and Run DMD
for ns = 1:length(splits)
    disp(['Running ns = ' num2str(ns)])
    freqs = freqs_unsplit;
    freqs(1:2:end-1) = freqs(1:2:end-1) - splits(ns)/2;
    freqs(2:2:end) = freqs(2:2:end) + splits(ns)/2;
    X = repmat(amps.',1,length(t)).*sin(freqs.' * t + phases.');

    all_om = zeros(nSteps*r,1);
    for k = 1:nSteps
        thisWind = (k-1)*stepSize + 1 : (k-1)*stepSize + wSize;
        xSample = X(:,thisWind);
        tSample = t(thisWind);
        tSample = tSample - tSample(1); %shift to start @ 0
        [Phi, Omega, b] = dmd(X,t,r);
        vs_res{ns,k}.Phi = Phi;
        vs_res{ns,k}.Omega = Omega;
        vs_res{ns,k}.b = b;
        vs_res{ns,k}.thisWind = thisWind;
        all_om((k-1)*r + 1 : k*r) = Omega;
    end
    vs_freqs{ns}.all_om = all_om;
    all_om_sq = all_om.*(conj(all_om));
    vs_freqs{ns}.all_om_sq = all_om_sq;

end
save('varied_split_res.mat','vs_freqs','vs_res','splits','r');
%% Plot Varied-Window Results
if exist('vs_freqs','var') == 0
    load('varied_split_res.mat');
end
for ns = 1:length(splits)
    all_om = vs_freqs{ns}.all_om;
    all_om_sq = vs_freqs{ns}.all_om_sq;
    abs_om = abs(all_om);
    figure('Position',[100 100 1100 600])
    subplot(1,2,1)
    ctrs = cell(2,1);
    ctrLim = max([max(abs(real(all_om))) max(abs(imag(all_om)))]);
    ctrs{1} = linspace(-ctrLim,ctrLim,length(freqs));
    ctrs{2} = ctrs{1};
    h3 = hist3([real(all_om) imag(all_om)],ctrs);
    h3_1 = h3';
    h3_1(size(h3,1) + 1, size(h3,2) + 1) = 0;
    xb = linspace(-100,100,size(h3,1)+1);
    yb = linspace(-100,100,size(h3,1)+1);
    hpc = pcolor(xb,yb,h3_1);
    hpc.ZData = ones(size(h3_1)) * -max(max(h3));
    colormap(hot)
    colorbar
    xlabel('Re[\omega]');
    ylabel('Im[\omega]');
    axis square
    title(['Freq. Split ' num2str(splits(ns))]);
    subplot(1,2,2)
    histogram(abs_om,length(freqs));
    yl = ylim;
    hold on
    plot([freqs; freqs],repmat(yl,length(freqs),1).','k','LineWidth',0.5);
    hold off
    xlim([min(abs_om), max(abs_om)]);
    xlabel('|\omega|');
    title(['Rank ' num2str(r) ' DMD']);
    
    outFile = ['./Figs/Varied_Splitting/Fig' num2str(ns) '.png'];
    saveas(gcf,outFile);
    
    
%     %% Plot Sample DMD Recon.
%     nTest = 1;
%     kTest = 5;
%     thisWind = vs_res{nTest,kTest}.thisWind;
%     Phi = vs_res{nTest,kTest}.Phi;
%     Omega = vs_res{nTest,kTest}.Omega;
%     b = vs_res{nTest,kTest}.b;
%     tTest = t(thisWind);
%     xTest = X(:,thisWind);
%     % xRecon = zeros(size(xTest));
%     time_dynamics = zeros(r,length(tTest));
%     for ts = 1:length(tTest)
%         time_dynamics(:,ts) = (b.*exp(Omega*tTest(ts)));
%     end
%     xRecon = Phi * time_dynamics;
%     figure
%     plot(tTest,xTest,'k',tTest,xRecon,'r')
end