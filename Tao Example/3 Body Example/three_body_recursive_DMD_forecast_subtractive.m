% Re-runs mwDMD leading up to forecast times as would need to be done in a
% realistic application. Recursive subtraction not yet implemented.


close all; clc
addpath(genpath('../optdmd-master'));

run('../optdmd-master/setup.m');


load('recursive_params.mat');
% if exist('pos','var') == 0
load('Three_Body_Data_Cartesian.mat');
% end

x = pos;
tspan_orig = tspan;
dt_orig = tspan(2)-tspan(1);

imode = 1;

n_test = 2;
tot_components = sum(rec_nComponents)  - sum(cell2mat(rec_composite_components));


dmd_hist_length_default = 10; % max # of historical DMD steps to use for forecasting
forecast_length_default = 5; % max # of window lengths to forecast


all_dt = zeros(n_recursions,1);
all_tBounds = zeros(n_recursions,2);
for j = 1:n_recursions
    nComponents = rec_nComponents(j);
    load(['mwDMD_sep_recon_recursion_' num2str(j,'%02.f') '.mat'], 'tspan');
    all_dt(j) = tspan(2)-tspan(1);
    all_tBounds(j,:) = [tspan(1) tspan(end)];
end
longest_window = max(rec_wSteps.*all_dt.');
longest_step = max(rec_stepSize.*all_dt.');

rng(111); %seed

all_test_times = [(max(all_tBounds(:,1))+longest_window+dmd_hist_length_default*longest_step) ...
    : min(all_dt) : (min(all_tBounds(:,2))-longest_window)];
all_test_times = all_test_times(randperm(length(all_test_times)));
all_test_times = sort(all_test_times(1:n_test));

all_test_steps = zeros(n_recursions,length(all_test_times));
for j = 1:n_recursions
    load(['mwDMD_sep_recon_recursion_' num2str(j,'%02.f') '.mat'],'tspan');
    for k = 1:length(all_test_times)
        [~,minStep] = min(abs(tspan - all_test_times(k)));
        all_test_steps(j,k) = minStep;
%         disp(['Pass ' num2str(j) ', Test ' num2str(k) ': Time ' num2str(test_times(k)) ' is at window ' num2str(minStep)])
    end
end

%% Establish Sampling Points For Each Recursion
burst_extents = [2 1 0 0]; % current step +/- this many
history_time = 10*max(rec_stepSize.*all_dt.');
sample_times = cell(n_test,n_recursions); %vectors of burst-sampling times
for j = 1:n_test
    test_time = all_test_times(j);
    query_times = cell(n_recursions,1); %vectors of all dense possible query times, of which the bursts will be a subset
    for q = 1:n_recursions
%         history_steps = 0:(floor(history_time/(rec_stepSize(q)*all_dt(q)))-1);
        query_times{q} = test_time + (-history_time:rec_stepSize(q)*all_dt(q):0);
    end
    sample_times{j,n_recursions} = query_times{n_recursions};
    for p = 1:(n_recursions-1)
        q = n_recursions-p;
        for k = 1:length(query_times{q+1})
            [~,closest_query] = min(abs(query_times{q+1}(k) - query_times{q}));
            sample_times{j,q} = [sample_times{j,q} ...
                query_times{q}(max([1 closest_query-burst_extents(q)]) : ...
                min([length(query_times{q}) closest_query+burst_extents(q)]))];
        end
    end
    
    % Plot query intervals for each pass
    figure
    for jj = 1:n_recursions
        plot(query_times{jj},jj*ones(size(query_times{jj})),'b.','MarkerSize',10)
        hold on
        plot(sample_times{j,jj},jj*ones(size(sample_times{j,jj})),'r.','MarkerSize',10)
        hold on
%         xlim([100000 103000])
        ylim([0 n_recursions]+0.5)
        ylabel('Recursive Iteration')
        set(gca,'YTick',1:n_recursions)
        xlabel('t')
    end
    
end

%% DMD Forecast
for j = 1:n_test
    fPlots(j) = figure;
end

for q = 1:n_recursions    
    load(['mwDMD_sep_recon_recursion_' num2str(q,'%02.f') '.mat']);
    % for "ground truth" comparison
    
    load(['mwDMD_params_recursion_' num2str(q,'%02.f') '.mat']);
    load(['km_centroids_recursion_' num2str(q,'%02.f') '.mat']);
    freq_meds = km_centroids.^(1/2);
    freq_meds = repelem(freq_meds,2); %duplicate each frequency for conjugate pair
    % global SVD
    [u,s,v] = svd(x,'econ');
    
    dt = all_dt(q);
    wSteps = rec_wSteps(q);
    stepSize = rec_stepSize(q);
    initialize_artificially = rec_initialize_artificially(q);
    r = rec_r(q);
    nComponents = rec_nComponents(q);
    tspan = all_tBounds(q,1) : all_dt(q) : all_tBounds(q,2);
    
    
    % FILTER DESIGN
    
    if q < n_recursions
        burst_interval = rec_wSteps(q+1) * all_dt(q+1);
    else
        burst_interval = rec_wSteps(q) * all_dt(q);
    end
    

    wind_t_rel = -2*dt*wSteps:dt:2*burst_interval;
%     wind_filt = exp(-(wind_t_rel-dt*wSteps/2).^2/(wSteps*dt/2)^2);
%     wind_filt = 0.05*tanh((wind_t_rel - dt*wSteps/2)/(dt*wSteps/4));
    mu = -dt*wSteps/2;
    sigma = dt*wSteps/2;
    tau = 4*dt*wSteps;
    wind_filt = (1 - normcdf((mu + (sigma.^2)/tau - wind_t_rel)./sigma)).*...
                   exp((mu + (sigma.^2)/(2*tau) - wind_t_rel)/tau)/tau;

    figure
    subplot(3,1,1)
    plot(wind_t_rel,wind_filt)
    hold on
    plot([0 0],ylim,'k--');
    hold on
    plot(wSteps*dt*[1 1],ylim,'k--')
    xlim([wind_t_rel(1) wind_t_rel(end)])
    title('Filter')
    all_wind_t = zeros(length(sample_times{j,q})*length(wind_t_rel),1);
    for jj = 1:length(sample_times{j,q})
        all_wind_t(((jj-1)*length(wind_t_rel)+1):jj*length(wind_t_rel)) = ...
            sample_times{j,q}(jj) + wind_t_rel;
    end
    
    all_wind_filt = repmat(wind_filt,1,length(sample_times{j,q})).';
    [t_join,filt_join] = sum_overlapping(all_wind_t,all_wind_filt);
    
    subplot(3,1,2)
    plot(t_join,filt_join);
    hold on
    for jj = 1:length(sample_times{j,q})
        plot(sample_times{j,q}(jj)*[1 1],ylim,'k--')
        hold on
    end
    title('Filter Sum')
    
    subplot(3,1,3)
    for jj = 1:length(sample_times{j,q})
        plot(all_wind_t(((jj-1)*length(wind_t_rel)+1):jj*length(wind_t_rel)), ...
            all_wind_filt(((jj-1)*length(wind_t_rel)+1):jj*length(wind_t_rel))./ ...
            interp1(t_join,filt_join,all_wind_t(((jj-1)*length(wind_t_rel)+1):jj*length(wind_t_rel))))
        hold on
    end
            
    % RECONSTRUCTION
    
        
    for j = 1:n_test
% for j = 3
        
%         dmd_hist_length = min([floor(all_test_steps(q,j)/rec_stepSize(q)), dmd_hist_length_default]);
        forecast_length = forecast_length_default*wSteps;
        if all_test_times(j) + forecast_length*dt > tspan(end)
            forecast_length = floor((tspan(end)-all_test_times(j))/(wSteps*dt)) * wSteps;
        end
        
%         t_test = all_test_times(j) + (-wSteps : forecast_length)*dt;
%         t_test_pad = all_test_times(j) + ((-wSteps-dmd_hist_length*stepSize):forecast_length)*dt;
        
        t_test = all_test_times(j) + (-history_time:dt:forecast_length*dt);
        t_test_pad = all_test_times(j) + ((-history_time-wSteps*dt):dt:forecast_length*dt);

        mr_res_local = cell(length(sample_times{j,q}),1);
        xr_sep_local = zeros(nVars,length(wind_t_rel),length(sample_times{j,q}),nComponents);
        t_local = zeros(length(wind_t_rel),length(sample_times{j,q}));
        [~,lower_step_bound] = min(abs(t_test_pad(1) - tspan_orig));
        [~,upper_step_bound] = min(abs(t_test_pad(end) - tspan_orig));
        steps_orig = lower_step_bound:upper_step_bound;
        clear('lower_step_bound','upper_step_bound');
        t_pad_orig = tspan_orig(steps_orig);
        x_local_orig = x(:,steps_orig);
        x_local = zeros(nVars,length(t_test_pad));
        for k = 1:nVars
            x_local(k,:) = interp1(t_pad_orig,x_local_orig(k,:),t_test_pad);
        end
        
        for k = 1:length(sample_times{j,q})
%             stepShift = (k-1)*stepSize;
%             windSteps = ((length(t_test_pad) - wSteps + 1) : length(t_test_pad)) - stepShift;
%             tShift = stepShift * dt;
            
            [~,lower_step_bound] = min(abs(t_test_pad - (sample_times{j,q}(k)-wSteps*dt)));
            [~,upper_step_bound] = min(abs(t_test_pad - (sample_times{j,q}(k))));
            windSteps = (lower_step_bound+1):upper_step_bound;
            clear('lower_step_bound','upper_step_bound');
            
            t_wind = t_test_pad(windSteps);
            x_wind = x_local(:,windSteps);
            
            % subtract mean
            c = mean(x_wind,2);
            x_wind = x_wind - repmat(c,1,size(x_wind,2));
            
            
            if (initialize_artificially == 1)
                % on first window, run uninitialized just to get phase angles
                if k == 1
                    [w, omega, b] = optdmd(x_wind,t_wind,r,imode,[],[],u);
                end
                [eSq, eInd] = sort(omega.*conj(omega)); %match order to that of freq_meds
                freq_angs = angle(omega(eInd));
                e_init = exp(sqrt(-1)*freq_angs) .* freq_meds;
                [w, omega, b] = optdmd(x_wind,t_wind,r,imode,[],e_init,u);
            else
                [w, omega, b] = optdmd(x_wind,t_wind,r,imode,[],[],u);
            end
            
            om_class = zeros(size(omega));
            om_sq = omega.*conj(omega);
            om_sq_dists  = (repmat(km_centroids.',r,1) - repmat(om_sq,1,nComponents)).^2;
            [~,om_class] = min(om_sq_dists,[],2);
            
            mr_res_local{k}.Omega = omega;
            mr_res_local{k}.w = w;
            mr_res_local{k}.b = b;
            mr_res_local{k}.c = c;
            mr_res_local{k}.t = t_wind;
            mr_res_local{k}.x = x_wind;
            mr_res_local{k}.om_class = om_class;
            
            if suppress_growth == 1
                omega(real(omega) > 0) = sqrt(-1) * imag(omega(real(omega) > 0));
            end
            
            this_recon_t = sample_times{j,q}(k) + wind_t_rel;
            t_local(:,k) = this_recon_t.';
            for nc = 1:nComponents
                xr_sep_local(:,:,k,nc) = w(:, om_class == nc)*diag(b(om_class == nc))*exp(omega(om_class == nc)*(this_recon_t-t_wind(1)));
                if nc == 1 % constant shift gets put in LF recon
                    xr_sep_local(:,:,k,nc) = xr_sep_local(:,:,k,nc) + c;
                end
                xr_sep_local(:,:,k,nc) = xr_sep_local(:,:,k,nc).*repmat(wind_filt,nVars,1);
%                 xr_sep{k,nc}(:,(k-1)*stepSize+1:(k-1)*stepSize+wSteps) = xr_sep{k,nc}(:,(k-1)*stepSize+1:(k-1)*stepSize+wSteps) + xr_sep_window{j};
            end
            
        
        end
        
        % Local Reconstruction
        xr_sep_local_join = cell(nComponents,1);
        
        for nc = 1:nComponents
            xr_sep_local_join{nc} = zeros(nVars,length(min(min(t_local)):dt:max(max(t_local))));
            for nv = 1:nVars
                [t_local_join,xr_sep_local_join{nc}(j,:)] = sum_overlapping(t_local,squeeze(xr_sep_local(nv,:,:,nc)));
                [~,filt_local_join] = sum_overlapping(t_local,repmat(wind_filt.',1,length(sample_times{j,q})));
            end
            xr_sep_local_join{nc} = xr_sep_local_join{nc}./repmat(filt_local_join.',nVars,1);
            plot(t_local_join,real(xr_sep_local_join{nc}))
            hold on
        end            
        
        

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
                xr_comp = w(:, om_class == k)*diag(b(om_class == k))*exp(omega(om_class == k)*(t_test-t_start));
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
                    omega = mr_res{histWind}.Omega;
                    if suppress_growth == 1
                        omega(real(omega) > 0) = sqrt(-1) * imag(omega(real(omega) > 0));
                    end
                    om_class = mr_res{histWind}.om_class;
                    t = mr_res{histWind}.t;
                    c = mr_res{histWind}.c;
                    t_start = mr_res{histWind}.t_start;
                    tShift = t-t(1);
                    
                    xr_hWind = w(:, om_class == k)*diag(b(om_class == k))* ...
                        exp(omega(om_class == k)*(t_hist-t_start));
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
