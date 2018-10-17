close all; clc

if exist('mr_res','var') == 0
    load('mwDMD_mr_res_i2.mat');
end
load('mwDMD_params.mat');
load('mwDMD_sep_recon.mat');
load('recon_gaussian_filter');

forecast_length = 5*wSteps;
nTestWind = 1; % number of randomized windows to test
rng(111); %seed
test_windows = randperm(nSlide-1000) + 500;
test_windows = sort(test_windows(1:nTestWind));

suppress_growth = 1; %eliminate all positive real components to frequencies
use_dmd_hist = 1; %dmd reconstruction using weighted contributions from last wSteps/stepSize steps

poly_fit_order = 5;

dmd_hist_length = floor(wSteps/stepSize);

dt = mr_res{1}.t(2) - mr_res{1}.t(1);

%% Collect Component-Wise Parameter Statistics
all_omega = cell(nComponents,1);
all_b = cell(nComponents,1);
sd_omega = cell(nComponents,1);
sd_b = cell(nComponents,1);
for k = 1:nSlide
    for j = 1:nComponents
        all_omega{j} = [all_omega{j}; mr_res{k}.Omega(mr_res{k}.om_class == j)];
        all_b{j} = [all_b{j}; mr_res{k}.Omega(mr_res{k}.om_class == j)];
    end
end
for j = 1:nComponents
    sd_omega{j} = [std(real(all_omega{j})) std(abs(imag(all_omega{j})))];
    sd_b{j} = [std(real(all_b{j})) std(abs(imag(all_b{j})))];
end

%% Forecast on test windows
for j = 1:length(test_windows)
    hist_windows = test_windows(j) + (-(dmd_hist_length-1):0);
    hist_windows = hist_windows(hist_windows >= 1); %in case of out-of-bounds indices
    nHist = length(hist_windows);
    
    t_test = mr_res{test_windows(j)}.t_start + dt*(0:forecast_length-1);
    steps_test = round(t_test/dt);
    
    % Build history matrices
    wHist = zeros(nHist,nVars,nVars);
    omegaHist = zeros(nHist,nVars);
    bHist = zeros(nHist,nVars);
    cHist = zeros(nHist,nVars);
    omClassHist = zeros(nHist,nVars);
    tStartHist = zeros(nHist,1);
    
    for ih = 1:nHist
        histWind = hist_windows(ih);
        wHist(ih,:,:) = mr_res{histWind}.w;
        omegaHist(ih,:) = mr_res{histWind}.Omega;
        bHist(ih,:) = mr_res{histWind}.b;
        cHist(ih,:) = mr_res{histWind}.c;
        omClassHist(ih,:) = mr_res{histWind}.om_class;
        tStartHist(ih) = mr_res{histWind}.t_start;
    end
    
    % Forecast
    xr = DMD_recon_forecast(t_test,wSteps,nComponents,recon_filter_sd,wHist,omegaHist,bHist,cHist,omClassHist,tStartHist,suppress_growth);
    figure
    for k = 1:nComponents
        xr_comp = xr{k};
        plot(t_test,xr_comp);
        hold on
    end
    
    % Plot
    figure
    for k = 1:nComponents
        subplot(nComponents,1,k)
        plot(tspan(steps_test),xr_sep{k}(:,steps_test),'k-')
        hold on
        plot(t_test,xr{k},'r-')
        hold on
        plot([t_test(wSteps) t_test(wSteps)],ylim,'k--')
        ylabel(['Recon: Component ' num2str(k)])
        xlim([t_test(1) t_test(end)]);
    end
end
