clear variables; close all; clc

nonlinear = 0; %toggle use of linear vs. nonlinear input data

if nonlinear == 0
    load('../raw_data_2_hiRes_linear.mat');
else
    load('../raw_data_2_hiRes.mat');
end

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
    if nonlinear == 0
        outFile = ['Data/varied_r_linear_' num2str(r)];
    else
        outFile = ['Data/varied_r_nonlinear_' num2str(r)];
    end
    [res_list, kmList, mr_res] = example2_tdDMD(x,TimeSpan,wSteps,nSplit,r,use_last_freq,thresh_pct,nComponents);    
    xr = example2_tdDMD_plot(x,TimeSpan,wSteps,nSplit,r,mr_res,kmList,thresh_pct,nComponents,fig_downsample,outFile);
    r_errs(j) = norm(x(:,1:fig_downsample:nSteps)-xr)/norm(x(:,1:fig_downsample:nSteps));
end
figure
plot(ranks,r_errs)
if nonlinear == 0
    saveas(gcf,'./Data/varied_rank_error_linear.fig');
else
    saveas(gcf,'./Data/varied_rank_error_nonlinear.fig');
end


%% Test Different Delays
if nonlinear == 0
    load('../raw_data_2_hiRes_linear.mat');
else
    load('../raw_data_2_hiRes.mat');
end
r = 16;
delays = floor(5 * 10.^(0:0.5:2.5));
delay_errs = zeros(size(delays));
for di = 1:length(delays)
    delaySteps = delays(di);
    xD = zeros(size(x,1)*nDelay, size(x,2)-(nDelay-1)*delaySteps);
    for j = 1:nDelay
        xD((j-1)*size(x,1)+1 : j*size(x,1),:) = x(: , (j-1)*delaySteps+1 : size(x,2)-(nDelay - j)*delaySteps);
    end
    if nonlinear == 0
        outFile = ['Data/varied_delay_linear_' num2str(di)];
    else
        outFile = ['Data/varied_delay_nonlinear' num2str(di)];
    end
    [res_list, kmList, mr_res] = example2_tdDMD(xD,TimeSpan,wSteps,nSplit,r,use_last_freq,thresh_pct,nComponents);    
    xr = example2_tdDMD_plot(xD,TimeSpan,wSteps,nSplit,r,mr_res,kmList,thresh_pct,nComponents,fig_downsample,outFile);
    delay_errs(di) = norm(xD(:,1:fig_downsample:nSteps)-xr)/norm(xD(:,1:fig_downsample:nSteps));
end
figure
semilogx(delays,delay_errs)
if nonlinear == 0
    saveas(gcf,'./Data/varied_delay_error_linear.fig');
else
    saveas(gcf,'./Data/varied_delay_error_nonlinear.fig');
end

%% Use x,y variables
if nonlinear == 0
    load('../raw_data_2_hiRes_linear.mat');
else
    load('../raw_data_2_hiRes.mat');
end
x = [x(1,:) + x(2,:);
    x(3,:).*cos(x(4,:));
    x(3,:).*sin(x(4,:));
    -x(1,:) + x(2,:).^3]; % transf. into x,y vars

r = 8;
x = x(3:4,:); %just y1, y2 (fast vars)

delays = floor(5 * 10.^(0:0.4:2.4));
delay_errs = zeros(size(delays));
for di = 1:length(delays)
    delaySteps = delays(di);
    xD = zeros(size(x,1)*nDelay, size(x,2)-(nDelay-1)*delaySteps);
    for j = 1:nDelay
        xD((j-1)*size(x,1)+1 : j*size(x,1),:) = x(: , (j-1)*delaySteps+1 : size(x,2)-(nDelay - j)*delaySteps);
    end
    if nonlinear == 0
        outFile = ['Data/varied_delay_xy_linear_' num2str(di)];
    else
        outFile = ['Data/varied_delay_xy_nonlinear' num2str(di)];
    end
    [res_list, kmList, mr_res] = example2_tdDMD(xD,TimeSpan,wSteps,nSplit,r,use_last_freq,thresh_pct,nComponents);    
    xr = example2_tdDMD_plot(xD,TimeSpan,wSteps,nSplit,r,mr_res,kmList,thresh_pct,nComponents,fig_downsample,outFile);
    title(['Delay = ' num2str((TimeSpan(2)-TimeSpan(1))*delaySteps)]);
    delay_errs(di) = norm(xD(:,1:fig_downsample:nSteps)-xr)/norm(xD(:,1:fig_downsample:nSteps));
end
figure
semilogx(delays,delay_errs)
if nonlinear == 0
    saveas(gcf,'./Data/varied_delay_xy_error_linear.fig');
else
    saveas(gcf,'./Data/varied_delay_xy_error_nonlinear.fig');
end
