clear variables; close all; clc

load('../raw_data_2_hiRes.mat');



% r = 12;
nComponents = 2;
wSteps = 11000;
nSplit = 20;
nSteps = wSteps * nSplit;
delaySteps = 200;
nDelay = 5;
use_last_freq = 0;
thresh_pct = 1.00;

fig_downsample = 10; %factor by which to downsample data for plotting
                     %should be an integer factor of wSteps

%% Add time delay embedding
xD = zeros(size(x,1)*nDelay, size(x,2)-(nDelay-1)*delaySteps);
for j = 1:nDelay
    xD((j-1)*size(x,1)+1 : j*size(x,1),:) = x(: , (j-1)*delaySteps+1 : size(x,2)-(nDelay - j)*delaySteps);
end

x = xD;
clear('xD');

%% Test Rank Dependence
ranks = 4:4:20;
r_errs = zeros(size(ranks));
for j = 1:length(ranks)
    r = ranks(j);
    outFile = ['Data/varied_r_' num2str(r)];
    [res_list, kmList, mr_res] = example2_tdDMD(x,TimeSpan,wSteps,nSplit,r,use_last_freq,thresh_pct,nComponents);    
    xr = example2_tdDMD_plot(x,TimeSpan,wSteps,nSplit,r,mr_res,kmList,thresh_pct,nComponents,fig_downsample,outFile);
    r_errs(j) = norm(x(:,1:fig_downsample:nSteps)-xr)/norm(x(:,1:fig_downsample:nSteps));
end
figure
plot(ranks,r_errs)
saveas(gcf,'./Data/varied_rank_error.fig')

%% Test Different Delays
load('../raw_data_2_hiRes.mat');
r = 16;
delays = floor(5 * 10.^(0:0.5:2.5));
delay_errs = zeros(size(ranks));
for di = 1:length(delays)
    delaySteps = delays(di);
    xD = zeros(size(x,1)*nDelay, size(x,2)-(nDelay-1)*delaySteps);
    for j = 1:nDelay
        xD((j-1)*size(x,1)+1 : j*size(x,1),:) = x(: , (j-1)*delaySteps+1 : size(x,2)-(nDelay - j)*delaySteps);
    end
    outFile = ['Data/varied_delay_' num2str(di)];
    [res_list, kmList, mr_res] = example2_tdDMD(xD,TimeSpan,wSteps,nSplit,r,use_last_freq,thresh_pct,nComponents);    
    xr = example2_tdDMD_plot(xD,TimeSpan,wSteps,nSplit,r,mr_res,kmList,thresh_pct,nComponents,fig_downsample,outFile);
    delay_errs(di) = norm(xD(:,1:fig_downsample:nSteps)-xr)/norm(xD(:,1:fig_downsample:nSteps));
end
figure
semilogx(delays,delay_errs)
saveas(gcf,'./Data/varied_delay_error.fig')