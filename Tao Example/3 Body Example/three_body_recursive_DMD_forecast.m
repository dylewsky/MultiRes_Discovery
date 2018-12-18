close all; clc
load('recursive_params.mat');
% if exist('all_mr_res','var') == 0
%     load('recursive_res.mat');
% end


n_test = 2;
use_dmd_hist = 1; %dmd reconstruction using weighted contributions from last wSteps/stepSize steps
reproject_b = 0;
ensemble_forecast = 1;
ensemble_alpha = 0.1;
% multipliers for parameters' standard deviations to use in random sampling
vary_params = [1 1 1]; % omega, w, b

nEnsemble = 20;
tot_components = sum(rec_nComponents)  - sum(cell2mat(rec_composite_components));

%% Establish test time domains
all_tBounds = zeros(n_recursions,2);
all_dt = zeros(n_recursions,1);
all_t_start = cell(n_recursions,1);
for j = 1:n_recursions
    nComponents = rec_nComponents(j);
    load(['mwDMD_sep_recon_recursion_' num2str(j,'%02.f') '.mat'], 'tspan');
    all_dt(j) = tspan(2)-tspan(1);
    all_tBounds(j,:) = [tspan(1) tspan(end)];
    
    load(['mwDMD_params_recursion_' num2str(j,'%02.f') '.mat'], 'initialize_artificially','nSlide');
    if initialize_artificially == 1
        load(['mwDMD_mr_res_recursion_' num2str(j,'%02.f') '_i2.mat'],'t_starts');
    else
        load(['mwDMD_mr_res_recursion_' num2str(j,'%02.f') '.mat'],'t_starts');
    end
    all_t_start{j} = t_starts;
end

rng(111); %seed

common_time = all_t_start{1};
for j = 2:n_recursions
    common_time = intersect(common_time,all_t_start{j});
end
if length(common_time) == 0
    disp('No common time found')
else
    common_time = sort(common_time(common_time > 0));
    common_time = common_time(1); %smallest common time
end
all_test_times = common_time:common_time:min(all_tBounds(:,2));
% all_test_times = max(all_tBounds(:,1)) : max(all_dt) : min(all_tBounds(:,2));
all_test_times = all_test_times(round(0.1*length(all_test_times)) : round(0.5*length(all_test_times)));
all_test_times = all_test_times(all_test_times > max(all_dt.*(rec_wSteps.')));
all_test_times = all_test_times(randperm(length(all_test_times)));
all_test_times = sort(all_test_times(1:n_test));

all_test_windows = zeros(n_recursions,length(all_test_times));
for j = 1:n_recursions
    t_starts = all_t_start{j};
    for k = 1:length(all_test_times)
        [~,minStep] = min(abs(all_t_start{j} - all_test_times(k)));
        all_test_windows(j,k) = minStep;
%         disp(['Pass ' num2str(j) ', Test ' num2str(k) ': Time ' num2str(test_times(k)) ' is at window ' num2str(minStep)])
    end
end

% shift all test times back by one window, so they line up at the end
% rather than the beginning
all_test_windows = all_test_windows - ones(size(all_test_windows));
all_test_times = all_test_times - repmat(all_dt.*(rec_wSteps.'),1,n_test);

%% Get Statistics for Ensemble Forecast (If Applicable)
all_stat_sds = cell(tot_components,1);
if ensemble_forecast == 1
    c_ind = 0;
    for q = 1:n_recursions
        load(['mwDMD_params_recursion_' num2str(q,'%02.f') '.mat']);
        if initialize_artificially == 1
            load(['mwDMD_mr_res_recursion_' num2str(q,'%02.f') '_i2.mat']);
        else
            load(['mwDMD_mr_res_recursion_' num2str(q,'%02.f') '.mat']);
        end
        for nc = 1:nComponents
            if ismember(nc,rec_composite_components{q})
                continue
            end
            c_ind = c_ind + 1;
            all_stat_sds{c_ind} = get_mr_res_stats(mr_res,nc);
        end
    end
end
save('mwDMD_mr_res_statistics.mat','all_stat_sds');
% load('mwDMD_mr_res_statistics.mat','all_stat_sds');

%% DMD Forecast
for j = 1:n_test
    fPlots(j) = figure;
end

forecast_lengths = [5 5 3 3 1 1];
dmd_hist_length_default = 10;

if ensemble_forecast == 1
    xr_ensemble = cell(nEnsemble,tot_components);
end

for q = 1:n_recursions
    load(['mwDMD_sep_recon_recursion_' num2str(q,'%02.f') '.mat']);
    load(['mwDMD_params_recursion_' num2str(q,'%02.f') '.mat']);
    if initialize_artificially == 1
        load(['mwDMD_mr_res_recursion_' num2str(q,'%02.f') '_i2.mat']);
    else
        load(['mwDMD_mr_res_recursion_' num2str(q,'%02.f') '.mat']);
    end
    test_windows = all_test_windows(q,:);
    

    dt = all_dt(q);

        
    for j = 1:length(test_windows)
% for j = 3
%         dmd_hist_length = min([floor(wSteps/stepSize), test_windows(j)]);
        dmd_hist_length = min([dmd_hist_length_default, test_windows(j)]);
        
%         wSteps = min_wSteps; %force all forecasting to be done on same time domain
        
        forecast_length = forecast_lengths(j)*wSteps;
        forecast_time = mr_res{test_windows(j)}.t_start + dt*(0:(forecast_length-1));
        
        [~,lower_step_bound] = min(abs(tspan - forecast_time(1)));
        [~,upper_step_bound] = min(abs(tspan - forecast_time(end)));
        steps_test = lower_step_bound:upper_step_bound;
        t_test = tspan(steps_test);
        t_test_pad = [t_test (t_test(end)+dt*(1:wSteps))];
        
        forecast_length = length(t_test);

%         t_test = mr_res{test_windows(j)}.t_start + dt*(0:forecast_length-1);
%         t_test_pad = mr_res{test_windows(j)}.t_start + dt*(0:forecast_length+wSteps-1);
%         if t_test(end) > tspan(end) % don't let it run past the end of time
%             t_test = mr_res{test_windows(j)}.t_start : dt : tspan(end);
%             t_test_pad = mr_res{test_windows(j)}.t_start : dt : (tspan(end) + dt*wSteps);
%             forecast_length = length(t_test);
%         end
% 
%         [C,steps_test,~] = intersect(tspan,t_test);


        w = mr_res{test_windows(j)}.w;
        b = mr_res{test_windows(j)}.b;
        Omega = mr_res{test_windows(j)}.Omega;

        if suppress_growth == 1
            Omega(real(Omega) > 0) = sqrt(-1) * imag(Omega(real(Omega) > 0));
        end
        om_class = mr_res{test_windows(j)}.om_class;
        t = mr_res{test_windows(j)}.t;
        c = mr_res{test_windows(j)}.c;
        t_start = mr_res{test_windows(j)}.t_start;
        tShift = t-t(1);

        %windows to use for long-memory reconstruction
        hist_windows = test_windows(j) + (-(dmd_hist_length-1):0);
        hist_windows = hist_windows(hist_windows >= 1); %in case of out-of-bounds indices
        nHist = length(hist_windows);

        figure(fPlots(j));


        for ne = 1:max([nEnsemble*ensemble_forecast, 1])
            % loop over ensemble iterations (if applicable)
            c_ind = 0;
            if q > 1
                for qi = 1:(q-1)
                    c_ind = c_ind + rec_nComponents(qi) - length(rec_composite_components{qi});
                end
            end
            w_mod = w;
            for k = 1:nComponents
                if ismember(k,rec_composite_components{q})
                    continue
                end
                c_ind = c_ind + 1;
                if ensemble_forecast == 1
                    mr_res_stats = all_stat_sds{c_ind};
                end


%                 vary_params = [1 1 1]; % omega, w, b
                if ensemble_forecast == 1
                    omega_ang = [abs(Omega), angle(Omega)];
                    omega_ang = omega_ang + ...
                        [mr_res_stats.omega_ang(1) * randn * vary_params(1), ...
                        mr_res_stats.omega_ang(2) * randn * vary_params(1)];
                    omega_mod = omega_ang(1) * exp(sqrt(-1)*omega_ang(2));
                    
                    b_mod = b + mr_res_stats.b(1) * randn(r,1) * vary_params(3) + ...
                        sqrt(-1) * mr_res_stats.b(2) * randn(r,1) * vary_params(3);
                    
                    
                    this_w = w(:, om_class == k);
                    if size(this_w,2) ~= 2
                        continue
                    end
                    this_w = this_w(:,1);
                    this_w_ang = [abs(this_w) angle(this_w)];
                    this_w_ang = this_w_ang + [mr_res_stats.w_ang(1) .* randn(nVars,1) * vary_params(2), ...
                        mr_res_stats.w_ang(2) .* randn(nVars,1) * vary_params(2)];
                    this_w = this_w_ang(:,1) .* exp(sqrt(-1) * this_w_ang(:,2));
                    this_w = [this_w conj(this_w)];
                    w_mod(:, om_class == k) = this_w;

                end
        
        
                subplot(tot_components,1,c_ind)

                % DMD Forecast
                if use_dmd_hist == 0
                    if ensemble_forecast == 0
                        xr_comp = w(:, om_class == k)*diag(b(om_class == k))*exp(Omega(om_class == k)*(t_test-t_start));
                    else
                        xr_comp = w_mod(:, om_class == k)*diag(b_mod(om_class == k))*exp(omega_mod(om_class == k)*(t_test-t_start));
                    end
                    if k == 1 % constant shift gets put in LF recon
                        xr_comp = xr_comp + c;
                    end
                    xr_comp = real(xr_comp);
                else
                    hWeights = (1/nHist) * ones(nHist,1);
    %                 hWeights = 1:nHist;

                    xn = zeros(1,(nHist-1)*stepSize+forecast_length);
                    xr_comp = zeros(nVars,length(xn));
                    t_hist = t_start + dt*( (-(nHist-1)*stepSize+1) : forecast_length);

                    for ih = 1:nHist
                        histWind = hist_windows(ih);
                        %                 hFilter = zeros(1,length(t_test));
                        hShift = (ih-1)*stepSize;

                        w = mr_res{histWind}.w;
                        b = mr_res{histWind}.b;


                        Omega = mr_res{histWind}.Omega;
                        if suppress_growth == 1
                            Omega(real(Omega) > 0) = sqrt(-1) * imag(Omega(real(Omega) > 0));
                        end
                        om_class = mr_res{histWind}.om_class;
                        t = mr_res{histWind}.t;
                        c = mr_res{histWind}.c;
                        t_start = mr_res{histWind}.t_start;
                        tShift = t-t(1);
                        
                        if ensemble_forecast == 1
                            omega_ang = [abs(Omega), angle(Omega)];
                            omega_ang = omega_ang + ...
                                [mr_res_stats.omega_ang(:,1) * randn(r,1) * vary_params(1), ...
                                mr_res_stats.omega_ang(:,2) * randn(r,1) * vary_params(1)];
                            omega_mod = omega_ang(:,1) .* exp(sqrt(-1)*omega_ang(:,2));

                            b_mod = b + mr_res_stats.b(1) * randn(r,1) * vary_params(3) + ...
                                sqrt(-1) * mr_res_stats.b(2) * randn(r,1) * vary_params(3);


                            this_w = w(:, om_class == k);
                            if size(this_w,2) ~= 2
                                continue
                            end
                            this_w = this_w(:,1);
                            this_w_ang = [abs(this_w) angle(this_w)];
                            this_w_ang = this_w_ang + [mr_res_stats.w_ang(1) .* randn(nVars,1) * vary_params(2), ...
                                mr_res_stats.w_ang(2) .* randn(nVars,1) * vary_params(2)];
                            this_w = this_w_ang(:,1) .* exp(sqrt(-1) * this_w_ang(:,2));
                            w_mod(:, om_class == k) = [this_w conj(this_w)];
%                             w_mod(:, om_class == k) = this_w;

                        end


                        if reproject_b == 1
                            reproj_steps = round(wSteps/4);
                            b_new = w(:, om_class == k) \ xr_sep{k}(:,steps_test((wSteps-reproj_steps):wSteps)) / exp(Omega(om_class == k)*(t_hist((wSteps-reproj_steps):wSteps)-t_start));
                            b_new = diag(b_new);
                            b(om_class == k) = b_new;
                        end


                        xr_hWind = w_mod(:, om_class == k)*diag(b_mod(om_class == k))* ...
                            exp(omega_mod(om_class == k)*(t_hist-t_start));
                        if k == 1 % constant shift gets put in LF recon
                            xr_hWind = xr_hWind + c;
                        end

                        xr_hWind = xr_hWind * hWeights(ih);
                        xn = xn + hWeights(ih);
    %                     xr_comp(:,((nHist-1)*stepSize - hShift + 1) : end) = ...
    %                         xr_comp(:,((nHist-1)*stepSize - hShift + 1) : end) + ...
    %                         xr_hWind(:,1:length(xr_comp(:,((nHist-1)*stepSize - hShift + 1) : end)));
                        xr_comp = xr_comp + xr_hWind;
    %                     xr_comp(:,wSteps-hShift+1 : end) = xr_comp(:,wSteps-hShift+1 : end) + ...
    %                         xr_hWind(:,1:end-wSteps+hShift);
    %                     xn(:,wSteps-hShift+1 : end) = xn(:,wSteps-hShift+1 : end) + hWeights(j);

                    end
                    xr_comp = xr_comp(:,(nHist-1)*stepSize + 1 : end); %cut off initial padding steps
                    xn = xn(:,(nHist-1)*stepSize + 1 : end);
                    xr_comp = real(xr_comp./repmat(xn,nVars,1)); %normalize
                end
                
                xr_ensemble{ne,c_ind} = xr_comp;

    %             comp_scale_loss = xr_comp./(xr_sep{k}(:,steps_test));
    %             comp_scale_loss = median(mean(comp_scale_loss));

                reconPlot = plot(t_test,xr_sep{k}(:,steps_test),'k-','DisplayName','Full DMD Recon.');
                if ensemble_forecast == 1
                    for np = 1:length(reconPlot)
                        reconPlot(np).Color(4) = ensemble_alpha;
                    end
                end
                ylabel(['Pass ' num2str(q) ', Component ' num2str(k)]);
                hold on
                plot(t_test,xr_comp,'r-','DisplayName','1-Window DMD Recon')
                hold on
                plot([t_test(wSteps) t_test(wSteps)],ylim,'k--','DisplayName','Boundary Between Window and Future')
                hold on
                xlim([t_test(1) t_test(end)])
                %         legend('Location','eastoutside');
            end
        end
    end
end
