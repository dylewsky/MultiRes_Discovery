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
all_omega_ang = cell(nComponents,1);
all_b = cell(nComponents,1);
all_w = cell(nComponents,1);
all_w_ang = cell(nComponents,2);
sd_omega = cell(nComponents,1);
sd_b = cell(nComponents,1);
sd_omega_ang = cell(nComponents,1);
sd_w_ang = cell(nComponents,1);
for k = 1:nSlide
    for j = 1:nComponents
        if nnz(mr_res{k}.om_class == j) ~= 2
            continue;
        end
        all_omega{j} = [all_omega{j}; mr_res{k}.Omega(mr_res{k}.om_class == j)];
        all_omega_ang{j} = [all_omega_ang{j}; ...
            [abs(mr_res{k}.Omega(mr_res{k}.om_class == j))...
            abs(angle(mr_res{k}.Omega(mr_res{k}.om_class == j)))]]; 
        all_b{j} = [all_b{j}; mr_res{k}.b(mr_res{k}.om_class == j)];
        this_w = mr_res{k}.w(:,mr_res{k}.om_class == j);
        all_w{j} = [all_w{j} this_w(:,1)]; %just keep one of conjugate pair
    end
end
for j = 1:nComponents
    all_w_ang{j,1} = abs(all_w{j});
    all_w_ang{j,2} = angle(all_w{j});
    
    sd_omega{j} = [std(real(all_omega{j})) std(abs(imag(all_omega{j})))];
    sd_omega_ang{j} = [std(all_omega_ang{j}(:,1)) std(all_omega_ang{j}(:,2))];
    sd_b{j} = [std(real(all_b{j})) std(abs(imag(all_b{j})))];
    sd_w_ang{j} = [std(all_w_ang{j,1},0,2) std(abs(all_w_ang{j,2}),0,2)];
end

% %% Forecast on test windows
% for j = 1:length(test_windows)
%     hist_windows = test_windows(j) + (-(dmd_hist_length-1):0);
%     hist_windows = hist_windows(hist_windows >= 1); %in case of out-of-bounds indices
%     nHist = length(hist_windows);
%     
%     t_test = mr_res{test_windows(j)}.t_start + dt*(0:forecast_length-1);
%     steps_test = round(t_test/dt);
%     
%     % Build history matrices
%     wHist = zeros(nHist,nVars,nVars);
%     omegaHist = zeros(nHist,nVars);
%     bHist = zeros(nHist,nVars);
%     cHist = zeros(nHist,nVars);
%     omClassHist = zeros(nHist,nVars);
%     tStartHist = zeros(nHist,1);
%     
%     for ih = 1:nHist
%         histWind = hist_windows(ih);
%         wHist(ih,:,:) = mr_res{histWind}.w;
%         omegaHist(ih,:) = mr_res{histWind}.Omega;
%         bHist(ih,:) = mr_res{histWind}.b;
%         cHist(ih,:) = mr_res{histWind}.c;
%         omClassHist(ih,:) = mr_res{histWind}.om_class;
%         tStartHist(ih) = mr_res{histWind}.t_start;
%     end
%     
%     % Forecast
%     xr = DMD_recon_forecast(t_test,wSteps,nComponents,recon_filter_sd,wHist,omegaHist,bHist,cHist,omClassHist,tStartHist,suppress_growth);
%     figure
%     for k = 1:nComponents
%         xr_comp = xr{k};
%         plot(t_test,xr_comp);
%         hold on
%     end
%     
%     % Plot
%     figure
%     for k = 1:nComponents
%         subplot(nComponents,1,k)
%         plot(tspan(steps_test),xr_sep{k}(:,steps_test),'k-')
%         hold on
%         plot(t_test,xr{k},'r-')
%         hold on
%         plot([t_test(wSteps) t_test(wSteps)],ylim,'k--')
%         ylabel(['Recon: Component ' num2str(k)])
%         xlim([t_test(1) t_test(end)]);
%     end
% end

%% Ensemble Forecast on Slow Scale
figure
% multipliers for parameters' standard deviations to use in random sampling
vary_params = [0 0 0]; % omega, w, b
% same, but for component biases
param_biases = [0 0 1]; % omega, w ,b
plot_component = 1; %which scale to plot
nEnsemble = 20;
ens_plots = cell(length(test_windows),nEnsemble);
ens_plot_alpha = 0.1;
for j = 1:length(test_windows)
    subplot(length(test_windows),1,j)
    
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
    
    plot(tspan(steps_test),xr_sep{plot_component}(:,steps_test),'k-','LineWidth',2)
    hold on
    
    % Apply stochastic variation
    xr_ensemble = cell(nEnsemble,nComponents);
    for q = 1:nEnsemble
        omegaHistMult = zeros(size(omegaHist));
        wHistMult = zeros(size(wHist));
        
        omegaMult_biases = zeros(nComponents,2); %r, theta
        wMult_biases = zeros(nVars,nVars,nComponents); %nVars nVars-element vectors for each component
        for p = 1:nComponents
            omegaMult_biases(p,:) = randn(1,2) .* sd_omega_ang{p} * param_biases(1);
            wMult_biases(:,:,p) = (param_biases(2) * randn(nVars,nVars) .* repmat(sd_w_ang{p}(:,1),1,nVars).') .* ...
                exp(sqrt(-1) * param_biases(2) * randn(nVars,nVars) .* repmat(sd_w_ang{p}(:,2),1,nVars).');
            if all(vecnorm(wMult_biases(:,:,p))) %if 0, don't normalize
                wMult_biases(:,:,p) = wMult_biases(:,:,p) ./ repmat(vecnorm(wMult_biases(:,:,p)),nVars,1);
            end
        end
        
%         disp(omegaMult_biases(1,:))
        
        for p = 1:nComponents
            omegaHistMult = omegaHistMult + ~(omClassHist-p) .* ...
                (   (ones(size(omegaHist)) + omegaMult_biases(p,1) + ...
                sd_omega_ang{p}(1) * vary_params(1) * randn(size(omegaHist))) .* ...
                exp(sqrt(-1) * (sd_omega_ang{p}(2) *  vary_params(1) * ...
                randn(size(omegaHist)) + repmat(omegaMult_biases(p,2),size(omegaHist)))) ...
                );
            
            wHistMult = wHistMult + permute(repmat(~(omClassHist-p),1,1,nVars),[1 3 2]) .* ...
                (   (ones(size(wHist)) + vary_params(2) * ...
                     permute(repmat(sd_w_ang{p}(:,1),1,nHist,nVars),[2 3 1]) .* ...
                     randn(size(wHist))) .* ...
                exp(sqrt(-1) *  vary_params(2) * permute(repmat(sd_w_ang{p}(:,2),1,nHist,nVars),[2 1 3]) .* ...
                randn(size(wHist))) ...
                );
        end
        
        for ih = 1:nHist
            for p = 1:nComponents
                if nnz(omClassHist(ih,:) == p) == 2 %check for conj pairs
                    omega_mult_pair = (ones(1,2) + sd_omega_ang{p}(1) * randn * vary_params(1) + omegaMult_biases(p,1))...
                        .* exp(sqrt(-1) * (sd_omega_ang{p}(2) * vary_params(1) * randn + omegaMult_biases(p,2)) * [1 -1]);
                    omegaHistMult(ih,omClassHist(ih,:) == p) = omega_mult_pair;
                    
                    w_mult_pair = (ones(nVars,2) +  repmat(sd_w_ang{p}(:,1).*randn(nVars,1),1,2) * vary_params(2))...
                        .* exp(sqrt(-1) * repmat(sd_w_ang{p}(:,2).*randn(nVars,1),1,2) * vary_params(2)); 
                    wHistMult(ih,:,omClassHist(ih,:) == p) = w_mult_pair;
                end
            end
        end
        
        % normalize and multiply by given coefficient
        wHistMult = vary_params(2) * wHistMult ./ repmat(vecnorm(wHistMult,2,2),1,nVars,1);
        for p = 1:nComponents
            % add (already normalized) bias, weighted by given coefficients
            wHistMult = wHistMult + param_biases(2) * ...
                permute(repmat(~(omClassHist-p),1,1,nVars),[1 3 2]) .* ...
                repmat(permute(wMult_biases(:,:,p),[3 1 2]),nHist,1,1);
        end
        % final normalization
        wHistMult = wHistMult ./ repmat(vecnorm(wHistMult,2,2),1,nVars,1);
        
        
        omegaHistMod = omegaHist .* omegaHistMult;
        
%         if nnz(isnan(wHistMod)) == 0
%             wHistMod = wHist .* wHistMult;
%             wHistMod = wHistMod ./ repmat(vecnorm(wHistMod,2,2),1,nVars,1);
%         else
%             disp('Using original w matrices')
%             wHistMod = wHist;
%         end
        wHistMod = wHist;

%         omegaHistMod = omegaHist + 0.2*randn*sqrt(-1);
        
        bHistMult = sin((q/nHist) * (1:nHist)).';
        bHistMult = bHistMult / sum(bHistMult);
        bHistMult = repmat(bHistMult,1,nVars);
        bHistMod = bHist .* bHistMult;
        
        xr_ensemble(q,:) = DMD_recon_forecast(t_test,wSteps,nComponents,recon_filter_sd,wHistMod,omegaHistMod,bHistMod,cHist,omClassHist,tStartHist,suppress_growth).';
        
        ens_plots{j,q} = plot(t_test,xr_ensemble{q,plot_component},'r-');
        for p = 1:nVars
            ens_plots{j,q}(p).Color(4) = ens_plot_alpha;
            hold on
        end
    end
    
    plot([t_test(wSteps) t_test(wSteps)],ylim,'k--')
    ylabel(['Recon: Component ' num2str(k)])
    xlim([t_test(1) t_test(end)]);
    
end
