clear variables; close all; clc

% testLabel = 'bimodal';
testLabel = 'pairwise';
% testLabel = 'uniform';

inData = [testLabel '_data.mat'];
load(inData);

if exist('fPeaks','var') == 0
    fPeaks = freqs;
end

%define which DMD function to use
dmd = @(x,t,r) reg_dmd(x,t,r);

%% Vary window size
windows = 10.^(1:0.5:0.5*floor(2*log10(length(t))));
target_nSteps = floor(length(t)/min(windows));

maxSteps = max(floor((length(t) - windows)./(floor((length(t) - windows)/target_nSteps))));

vw_res = cell(length(windows),maxSteps);
vw_freqs = cell(length(windows),1);

r = N;

for n = 1:length(windows)
    disp(['Running n = ' num2str(n)])
    wSize = windows(n);
    stepSize = 10 * floor((length(t) - wSize)/target_nSteps); %~10x larger than smallest window
    nSteps = floor((length(t) - wSize)/stepSize);
    all_om = zeros(nSteps*r,1);
    for k = 1:nSteps
        thisWind = (k-1)*stepSize + 1 : (k-1)*stepSize + wSize;
        xSample = X(:,thisWind);
        tSample = t(thisWind);
        tSample = tSample - tSample(1); %shift to start @ 0
        [Phi, Omega, b] = dmd(X,t,r);
        vw_res{n,k}.Phi = Phi;
        vw_res{n,k}.Omega = Omega;
        vw_res{n,k}.b = b;
        vw_res{n,k}.thisWind = thisWind;
        all_om((k-1)*r + 1 : k*r) = Omega;
    end
    vw_freqs{n}.all_om = all_om;
    all_om_sq = all_om.*(conj(all_om));
    vw_freqs{n}.all_om_sq = all_om_sq;
end
save([testLabel '_varied_window_res.mat'],'vw_freqs','vw_res','windows','r');

%% Plot Varied-Window Results
if exist('vw_freqs','var') == 0
    load([testLabel '_varied_window_res.mat']);
end
r = N;
for n = 1:length(windows)
    all_om = vw_freqs{n}.all_om;
    all_om_sq = vw_freqs{n}.all_om_sq;
    abs_om = abs(all_om);
    figure('Position',[100 50 1100 600])
    subplot(1,2,1)
    ctrs = cell(2,1);
    ctrLim = max([max(abs(real(all_om))) max(abs(imag(all_om)))]);
    ctrs{1} = linspace(-ctrLim,ctrLim,N);
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
    title([num2str(windows(n)) '-Step Window']);
    subplot(1,2,2)
    histogram(abs_om,N);
    yl = ylim;
    hold on
    plot([fPeaks; fPeaks],repmat(yl,length(fPeaks),1).','k');
    hold off
    xlabel('|\omega|');
    title(['Rank ' num2str(r) ' DMD']);
    outFile = ['./Figs/Varied_Window/' testLabel '_Fig' num2str(n) '.png'];
    saveas(gcf,outFile)
end

%% Plot Sample DMD Recon.
nTest = 1;
kTest = 5;
thisWind = vw_res{nTest,kTest}.thisWind;
Phi = vw_res{nTest,kTest}.Phi;
Omega = vw_res{nTest,kTest}.Omega;
b = vw_res{nTest,kTest}.b;
tTest = t(thisWind);
xTest = X(:,thisWind);
% xRecon = zeros(size(xTest));
time_dynamics = zeros(r,length(tTest));
for ts = 1:length(tTest)
    time_dynamics(:,ts) = (b.*exp(Omega*tTest(ts)));
end
xRecon = Phi * time_dynamics;
figure
plot(tTest,xTest,'k',tTest,xRecon,'r')

%% Visualize Windows
figure('Position',[100 50 900 600])
plot(t,X([1 33 67 100],:))
hold on
yMin = min(min(X([1 33 67 100],:))) - 1;
for n = 1:length(windows)
    wSize = floor(windows(n));
    yPos = yMin + 0.1 * n;
    plot([t(1) t(wSize)],[yPos yPos],'r','LineWidth',3);
    hold on
end
title('Input Signal & Window Sizes');
xlim([0 t(floor(max(windows))) + 1])
saveas(gcf,['./Figs/Varied_Window/' testLabel '_wSizes.png']);

%% Initialize Varied Rank
wSize = 200;
stepSize = 10;
nSteps = floor((length(t) - wSize)/stepSize);

ranks = [100 90 80 70 50 25 10];

%% Varied-Rank DMD

vr_res = cell(length(ranks),nSteps);
vr_freqs = cell(length(ranks),1);

for n = 1:length(ranks)
    disp(['Running n = ' num2str(n)])
    r = ranks(n);
    all_om = zeros(nSteps*r,1);
    for k = 1:nSteps
        thisWind = (k-1)*stepSize + 1 : (k-1)*stepSize + wSize;
        xSample = X(:,thisWind);
        tSample = t(thisWind);
        tSample = tSample - tSample(1); %shift to start @ 0
        [Phi, Omega, b] = dmd(X,t,r);
        vr_res{n,k}.Phi = Phi;
        vr_res{n,k}.Omega = Omega;
        vr_res{n,k}.b = b;
        vr_res{n,k}.thisWind = thisWind;
        all_om((k-1)*r + 1 : k*r) = Omega;
    end
    vr_freqs{n}.all_om = all_om;
    all_om_sq = all_om.*(conj(all_om));
    vr_freqs{n}.all_om_sq = all_om_sq;
end
save([testLabel '_varied_rank_res.mat'],'vr_freqs','vr_res','wSize','stepSize','nSteps','ranks');

%% Plot Varied-Rank Results
if exist('vr_freqs','var') == 0
    load([testLabel '_varied_rank_res.mat']);
end
for n = 1:length(ranks)
    all_om = vr_freqs{n}.all_om;
    all_om_sq = vr_freqs{n}.all_om_sq;
    abs_om = abs(all_om);
    figure('Position',[100 50 1100 600])
    subplot(1,2,1)
    ctrs = cell(2,1);
    ctrLim = max([max(abs(real(all_om))) max(abs(imag(all_om)))]);
    ctrs{1} = linspace(-ctrLim,ctrLim,N);
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
    title([num2str(wSize) '-Step Window']);
    subplot(1,2,2)
    histogram(abs_om,N);
    yl = ylim;
    hold on
    plot([fPeaks; fPeaks],repmat(yl,length(fPeaks),1).','k');
    hold off
    xlabel('|\omega|');
    title(['Rank ' num2str(ranks(n)) ' DMD']);
    outFile = ['./Figs/Varied_Rank/' testLabel '_Fig' num2str(n) '.png'];
    saveas(gcf,outFile)
end

%% Plot Varied-Rank DMD Recon.
kTest = 5;
wSize = 200;
figure('Position',[100 50 1100 600])
x_subset = [1 33 67 100];
for nTest = 1:length(ranks)
    r = ranks(nTest);
    thisWind = vr_res{nTest,kTest}.thisWind;
    Phi = vr_res{nTest,kTest}.Phi;
    Omega = vr_res{nTest,kTest}.Omega;
    b = vr_res{nTest,kTest}.b;
    tTest = t(thisWind);
    xTest = X(x_subset,thisWind);
    % xRecon = zeros(size(xTest));
    time_dynamics = zeros(r,length(tTest));
    for ts = 1:length(tTest)
        time_dynamics(:,ts) = (b.*exp(Omega*tTest(ts)));
    end
    xRecon = Phi * time_dynamics;
    subplot(3,3,nTest)
    plot(tTest,xTest,'k',tTest,xRecon(x_subset,:),'r','LineWidth',2)
    title(['r = ' num2str(r)])
    xlim([tTest(1) tTest(end)])
end
saveas(gcf,['./Figs/Varied_Rank/' testLabel '_Recons.png']);
