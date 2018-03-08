clear variables; close all; clc

testLabel = 'pairwise';
% testLabel = 'uniform';

inData = [testLabel '_data.mat'];
load(inData);

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
for n = 1:length(windows)
    all_om = vw_freqs{n}.all_om;
    all_om_sq = vw_freqs{n}.all_om_sq;
    figure
    subplot(1,2,1)
    ctrs = cell(2,1);
    ctrLim = max([max(abs(real(all_om))) max(abs(imag(all_om)))]);
    ctrs{1} = linspace(-ctrLim,ctrLim,length(freqs));
    ctrs{2} = ctrs{1};
    hist3([real(all_om) imag(all_om)],ctrs)
    xlabel('Re[\omega]');
    ylabel('Im[\omega]');
    axis square
    title([num2str(windows(n)) '-Step Window']);
    subplot(1,2,2)
    histogram(all_om_sq,length(freqs));
    yl = ylim;
    hold on
    plot([freqs.^2; freqs.^2],repmat(yl,length(freqs),1).','k');
    hold off
    xlabel('|\omega|^2');
    title(['Rank ' num2str(r) ' DMD']);
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

