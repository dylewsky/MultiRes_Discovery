close all; clc
load('recursive_params.mat');
% if exist('all_mr_res','var') == 0
%     load('recursive_res.mat');
% end


n_test = 2;
use_dmd_hist = 0; %dmd reconstruction using weighted contributions from last wSteps/stepSize steps

    

tot_components = sum(rec_nComponents)  - sum(cell2mat(rec_composite_components));

% all_t_start = cell(n_recursions,1);
% for j = 1:n_recursions
%     all_t_start{j} = all_mr_res{j}{1}.t_start : (all_mr_res{j}{2}.t_start - ...
%         all_mr_res{j}{1}.t_start) : all_mr_res{j}{end}.t_start;
% end

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


%% DMD Forecast
for j = 1:n_test
    fPlots(j) = figure;
end

forecast_lengths = [5 5 3 3 1.5 1.5];

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
        dmd_hist_length = min([floor(wSteps/stepSize), test_windows(j)]);
        
%         wSteps = min_wSteps; %force all forecasting to be done on same time domain
        
        forecast_length = forecast_lengths*(j)*wSteps;


        t_test = mr_res{test_windows(j)}.t_start + dt*(0:forecast_length-1);
        t_test_pad = mr_res{test_windows(j)}.t_start + dt*(0:forecast_length+wSteps-1);
        if t_test(end) > mr_res{end}.t(end) % don't let it run past the end of time
            t_test = mr_res{test_windows(j)}.t_start : dt : mr_res{end}.t(end);
            t_test_pad = mr_res{test_windows(j)}.t_start : dt : (mr_res{end}.t(end) + dt*wSteps);
            forecast_length = length(t_test);
        end

        steps_test = round(t_test/dt);

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
        c_ind = 0;
        if q > 1
            for qi = 1:(q-1)
                c_ind = c_ind + rec_nComponents(qi) - length(rec_composite_components{qi});
            end
        end
        for k = 1:nComponents
            if ismember(k,rec_composite_components{q})
                continue
            end
            c_ind = c_ind + 1;
                
            subplot(tot_components,1,c_ind)
            
            % DMD Forecast
            if use_dmd_hist == 0
                xr_comp = w(:, om_class == k)*diag(b(om_class == k))*exp(Omega(om_class == k)*(t_test-t_start));
                if k == 1 % constant shift gets put in LF recon
                    xr_comp = xr_comp + c;
                end
                xr_comp = real(xr_comp);
            else
                hWeights = (1/nHist) * ones(nHist,1);
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
                    
                    xr_hWind = w(:, om_class == k)*diag(b(om_class == k))* ...
                        exp(Omega(om_class == k)*(t_hist-t_start));
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
            
%             comp_scale_loss = xr_comp./(xr_sep{k}(:,steps_test));
%             comp_scale_loss = median(mean(comp_scale_loss));
            
            plot(t_test,xr_sep{k}(:,steps_test),'k-','DisplayName','Full DMD Recon.')
            title(['DMD Forecast at t = ' num2str(t_test(1)) ': Pass ' num2str(q) ', Component ' num2str(k)]);
            hold on
            plot(t_test,xr_comp,'r-','DisplayName','1-Window DMD Recon')
            hold on
            plot([t_test(wSteps) t_test(wSteps)],ylim,'k--','DisplayName','Boundary Between Window and Future')
            hold off
            xlim([t_test(1) t_test(end)])
            %         legend('Location','eastoutside');
        end
    end
end
