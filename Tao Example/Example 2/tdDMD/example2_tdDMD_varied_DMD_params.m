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

%% Add time delay embedding
xD = zeros(size(x,1)*nDelay, size(x,2)-(nDelay-1)*delaySteps);
for j = 1:nDelay
    xD((j-1)*size(x,1)+1 : j*size(x,1),:) = x(: , (j-1)*delaySteps+1 : size(x,2)-(nDelay - j)*delaySteps);
end

x = xD;
clear('xD');

%% Test Rank Dependence
ranks = 4:4:20;
errs = zeros(size(ranks));
for j = 1:length(ranks)
    r = ranks(j);
    outFile = ['Data/varied_r_' num2str(r)];
    [res_list, kmList, mr_res] = example2_tdDMD(x,TimeSpan,wSteps,nSplit,r,use_last_freq,thresh_pct,nComponents);
    xr = example2_tdDMD_plot(x,TimeSpan,wSteps,nSplit,r,mr_res,kmList,thresh_pct,nComponents,outFile);
    errs(j) = norm(x(:,1:nSteps)-xr)/norm(x(:,1:nSteps));
end
figure
plot(ranks,errs)