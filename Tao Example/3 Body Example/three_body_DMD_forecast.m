close all; clc

if exist('mr_res','var') == 0
    load('mwDMD_mr_res_i2.mat');
end
load('mwDMD_params.mat');
load('mwDMD_sep_recon.mat');
load('recon_gaussian_filter');

forecast_length = 5*wSteps;
rng(111); %seed
test_windows = randperm(nSlide-1000) + 500;
test_windows = sort(test_windows(1:5));

suppress_growth = 1; %eliminate all positive real components to frequencies
use_dmd_hist = 1; %dmd reconstruction using weighted contributions from last wSteps/stepSize steps

poly_fit_order = 5;

dmd_hist_length = floor(wSteps/stepSize);

dt = mr_res{1}.t(2) - mr_res{1}.t(1);
%% DMD vs. Polynomial vs HAVOK Forecast
for j = 1:length(test_windows)
% for j = 3
    t_test = mr_res{test_windows(j)}.t_start + dt*(0:forecast_length-1);
    t_test_pad = mr_res{test_windows(j)}.t_start + dt*(0:forecast_length+wSteps-1);
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
    
    fPlots(j) = figure;
    for k = 1:nComponents
        
        % DMD Forecast
        if use_dmd_hist == 0
            xr_comp = w(:, om_class == k)*diag(b(om_class == k))*exp(Omega(om_class == k)*(t_test-t_start));
            if k == 1 % constant shift gets put in LF recon
                xr_comp = xr_comp + c; 
            end
            xr_comp = real(xr_comp);
        else
            xr_comp = zeros(nVars,length(t_test)+wSteps);
            %pad with extra wSteps to allow reconstructions to be backshifted accordingly
            xn = zeros(1,length(t_test));
            
            % smallest value gaussian can take within window
            % (i.e. min of tail of oldest history window which still has
            % some overlap with the current window)
            filter_wSteps_min = exp(-(1.5*wSteps).^2/recon_filter_sd^2);
            
            % Build generic filter
            hFilter = exp(-((-stepSize + 1 : (length(t_test)+wSteps)) - (wSteps/2)).^2/recon_filter_sd^2);
%             hFilter = vpa(exp(-((-stepSize + 1 : (length(t_test)+wSteps)) - (wSteps/2)).^2/recon_filter_sd^2));

            hFilterLog = -((-stepSize + 1 : (length(t_test)+wSteps)) - (wSteps/2)).^2/recon_filter_sd^2;
%             % step at which gaussian is replaced by exponential:
% %                 exp_tau = 500; % # steps for tail to decay by 1/e
% %                 cross_step = round(recon_filter_sd^2 / (2 * exp_tau) + wSteps/2);
%             cross_step = stepSize + wSteps + recon_filter_sd ;
%             exp_tau = recon_filter_sd^2 / (2 * (cross_step -  (wSteps/2)));
% 
%             if cross_step <= wSteps/2 
%                 disp('Error: Tau too large')
%             elseif cross_step >= length(t_test)
%                 disp('Error: Tau too small')
%             end
%             exp_coeff = exp(recon_filter_sd^2/(4*exp_tau) + (wSteps/2)/exp_tau);
%             filter_asymptote = 0.001 * hFilter(cross_step);
            % must be small to maintain continuity of val + deriv
            % at cross_step
%             hFilter(cross_step : end) = exp_coeff * exp(-(1/exp_tau)*(cross_step : length(hFilter))) + filter_asymptote;
%             hFilter(cross_step : end) = (hFilter(cross_step) - filter_asymptote) * exp(-(1/exp_tau)*(0 : length(hFilter) - cross_step)) + filter_asymptote;
            
%             % Eliminate NaNs from filter
%             filt_min = 10^(-20);
%             num_iter = 0;
%             while nnz(isnan(hFilter)) ~= 0
%                 num_iter = num_iter + 1;
%                 for fj = 1:length(hFilter)
%                     if hFilter(fj) <= filt_min
%                         if fj < cross_step
%                             disp('Error: NaN before cross_step')
%                         else
%                             hFilter(fj : end) = (hFilter(fj) - filter_asymptote) * exp(-(1/exp_tau)*(0 : length(hFilter) - fj)) + filter_asymptote;
%                         end
%                         break
%                     end
%                 end
%                 if num_iter > 100
%                     disp('Failed to purge NaNs')
%                     break
%                 end
%             end
%             clear('num_iter','filt_min');

            allFilters = zeros(nHist+1,length(t_test));
%             allFiltersLog = zeros(nHist+1,length(t_test));
            for ih = 1:nHist+1 % 1 extra element that will dominate at large t
                hShift = (ih-1)*stepSize; %extra element gets negative shift
                allFilters(ih,:) = hFilter( hShift + 1 : length(hFilter) - wSteps - stepSize + hShift);
%                 allFiltersLog(ih,:) = hFilterLog( hShift + 1 : length(hFilterLog) - wSteps - stepSize + hShift);
                if k == 1 && j == 1 && mod(ih,10) == 1
                    figure
                    plot(t_test,allFilters(ih,:))
                    title(['Filter for ih = ' num2str(ih)])
%                     hold on
%                     plot([t_test(cross_step) t_test(cross_step)],[0 1],'r--','DisplayName','Piecewise Junction of Filter')
                    hold off
                    set(0, 'currentfigure', fPlots(j)); %change back to main figure
                end
            end
            pct_contributions = allFilters./(repmat(sum(allFilters),nHist+1,1));
%             pct_contributions_rel = zeros(size(allFiltersLog));
%             pct_contributions_rel(1,:) = 1; %reference value
%             u_ref = wSteps/2; %gaussian mean for reference filter
%             for ih = 2:nHist+1
%                 u_ih = wSteps/2 - (ih - 1);
%                 pct_contributions_rel(ih,:) = exp((u_ref^2 - u_ih^2)/recon_filter_sd^2) * ...
%                     exp(-2*t_test*(u_ref - u_ih)/recon_filter_sd^2);
%             end
%             pct_contributions_rel = pct_contributions_rel(2:end,:);
%             pct_contributions_rel = pct_contributions_rel./repmat(vecnorm(pct_contributions_rel),nHist,1);
%             figure
%             plot(t_test,pct_contributions_rel);

            % cut off extra element
            allFilters = allFilters(2:nHist+1,:);
            pct_contributions = pct_contributions(2:nHist+1,:);
            % add DC asymptote contribution
            pct_contributions = pct_contributions + 1/nHist;
            % renormalize
            pct_contributions = pct_contributions./repmat(sum(pct_contributions,1),nHist,1);
            
            
            % Interpolate over NaNs
            nan_cutoff = find(isnan(pct_contributions(1,:)));
            nan_cutoff = nan_cutoff(1);
            nan_cutoff = nan_cutoff - floor(wSteps/2); %cut off extra to omit rounding error noise near end
            pct_contributions(:,nan_cutoff:end) = NaN;
            
            
            pct_contributions(:,end) = 1/nHist *ones(nHist,1);
            pct_contributions = fillmissing(pct_contributions.','spline').';
            
            %normalize again
            pct_contributions = pct_contributions./repmat(sum(pct_contributions,1),nHist,1);
            
            allFilters = pct_contributions;
            % just use fractional contributions directly as conv. filter
            
            % Create reconstruction
            
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
                
                xr_hWind = w(:, om_class == k)*diag(b(om_class == k))*exp(Omega(om_class == k)*(t_test_pad-t_start));
                if k == 1 % constant shift gets put in LF recon
                    xr_hWind = xr_hWind + c; 
                end
                
                thisFilter = allFilters(ih,:);
                
                xr_comp(:,wSteps-hShift+1 : end) = xr_comp(:,wSteps-hShift+1 : end) + xr_hWind(:,1:end-wSteps+hShift);
                xr_comp(:,wSteps+1 : end) = xr_comp(:,wSteps+1 : end).*repmat(thisFilter,nVars,1);
                %no need to apply filter to first wSteps, those will be cut
                %off anyway
                xn = xn + thisFilter;
                
            end
            xr_comp = xr_comp(:,wSteps+1 : end); %cut off iniitial padding steps
            xr_comp = real(xr_comp./repmat(xn,nVars,1)); %normalize
            if k == 1 && j == 1 
                figure
                plot(t_test,pct_contributions)
                title('Window contribution fractions over time')
                xlim([t_test(1) t_test(end)])
                hold on
                plot([t_test(wSteps) t_test(wSteps)],[0 1],'k--');
%                 hold on
%                 plot([t_test(cross_step) t_test(cross_step)],[0 1],'r--','DisplayName','Piecewise Junction of Filter')
                hold off
                set(0, 'currentfigure', fPlots(j)); %change back to main figure
            end
        end
        subplot(nComponents,4,4*(k-1)+1)
        plot(t_test,xr_sep{k}(:,steps_test),'k-','DisplayName','Full DMD Recon.')
        title(['DMD Forecast at t = ' num2str(t_test(1)) ': Component ' num2str(k)]);
        hold on
        plot(t_test,xr_comp,'r-','DisplayName','1-Window DMD Recon')
        hold on
        plot([t_test(wSteps) t_test(wSteps)],ylim,'k--','DisplayName','Boundary Between Window and Future')
        hold off
        xlim([t_test(1) t_test(end)])
%         legend('Location','eastoutside');

        % Polynomial Forecast
        
        t_train = t_test(1:wSteps);
        xr_train = xr_sep{k}(:,steps_test(1:wSteps));
        xr_mean = mean(xr_train,2);
        xr_train = xr_train - xr_mean;
        pc = cell(nVars,1); %polynomial fit coefficients for each variable
        poly_forecast = cell(nVars,1);
        for vi = 1:nVars
            pc{vi} = polyfit(t_train,xr_train(vi,:),poly_fit_order);
            poly_forecast{vi} = polyval(pc{vi},t_test) + xr_mean(vi);
        end
        yBounds = ylim;
        subplot(nComponents,4,4*(k-1)+2)
        plot(t_test,xr_sep{k}(:,steps_test),'k-','DisplayName','Full DMD Recon.')
        title(['Polynomial Forecast at t = ' num2str(t_test(1)) ': Component ' num2str(k)]);
        hold on
        for vi = 1:nVars
            plot(t_test,poly_forecast{vi},'r-','DisplayName','1-Window Poly Recon')
            hold on
        end
        plot([t_test(wSteps) t_test(wSteps)],ylim,'k--','DisplayName','Boundary Between Window and Future')
        hold off
        xlim([t_test(1) t_test(end)])
        ylim(yBounds)
        
        %Global HAVOK Forecast
        
        load(['MR_HAVOK_model_' num2str(k) '.mat']);
        x0 = zeros(nDelay*size(xr_sep{k},1),1);
%         t0 = zeros(nDelay,1);
        for nd = 1:nDelay
            delayInds = ((nd-1)*size(xr_sep{k},1) + 1) : nd*size(xr_sep{k},1);
            x0(delayInds) = xr_sep{k}(:,steps_test(wSteps) - (nd-1)*delaySteps);
%             t0(nd) = tspan(steps_test(wSteps) - (nd-1)*delaySteps);
        end
        v0 = S^(-1) * U.' * x0; %transform into eigenbasis used by HAVOK model
        v0 = v0(1:r-1); %...up to rank used by HAVOK model

        [tHAV, vHAV] = ode45(@(tHAV,vHAV)HAVOK_rhs(vHAV,A),t_test(wSteps-1:end),v0);
        xHAV = U(:,1:r-1) * S(1:r-1,1:r-1) * vHAV.';
        xHAV = xHAV(1:nVars,:);
        subplot(nComponents,4,4*(k-1)+3)
        plot(t_test,xr_sep{k}(:,steps_test),'k-','DisplayName','Full DMD Recon.')
        title(['Global HAVOK Forecast at t = ' num2str(t_test(1)) ': Component ' num2str(k)]);
        hold on
        plot(tHAV,xHAV,'r-','DisplayName','HAVOK Forecast')
        hold on
        plot([t_test(wSteps) t_test(wSteps)],ylim,'k--','DisplayName','Boundary Between Window and Future')
        hold off
        xlim([t_test(1) t_test(end)])
        ylim(yBounds)
        
        % Local HAVOK Forecast
        
        % EIGEN-TIME DELAY COORDINATES
        Hloc = zeros(nVars*nDelay,length(steps_test)-(nDelay-1)*delaySteps);
        for hk=1:nDelay
            delayInds = ((hk-1)*nVars + 1) : hk*nVars;
            Hloc(delayInds,:) = xr_sep{k}(:,(hk-1)*delaySteps+1:(length(steps_test)-(nDelay)*delaySteps + hk*delaySteps));
        end

        [Uloc,Sloc,Vloc] = svd(Hloc,'econ'); % Eigen delay coordinates

        % COMPUTE DERIVATIVES (4TH ORDER CENTRAL DIFFERENCE)
        dVloc = zeros(length(Vloc)-5,r);
        for i=3:length(Vloc)-3
            for hk=1:r
                dVloc(i-2,hk) = (1/(12 * dt)) * (-Vloc(i+2,hk)+8 * Vloc(i+1,hk)-8 * Vloc(i-1,hk)+Vloc(i-2,hk));
            end
        end
        % trim first and last two that are lost in derivative
        Vloc = Vloc(3:end-3,1:r);

        % BUILD HAVOK REGRESSION MODEL ON TIME DELAY COORDINATES
        XiLoc = Vloc\dVloc;
        Aloc = XiLoc(1:r-1,1:r-1)';
        Bloc = XiLoc(end,1:r-1)';
        
        v0loc = Sloc^(-1) * Uloc.' * x0; %transform into eigenbasis used by HAVOK model
        v0loc = v0loc(1:r-1); %...up to rank used by HAVOK model

        [tHAVloc, vHAVloc] = ode45(@(tHAVloc,vHAVloc)HAVOK_rhs(vHAVloc,Aloc),t_test(wSteps-1:end),v0loc);
        xHAVloc = Uloc(:,1:r-1) * Sloc(1:r-1,1:r-1) * vHAVloc.';
        xHAVloc = xHAVloc(1:nVars,:);
        subplot(nComponents,4,4*(k-1)+4)
        plot(t_test,xr_sep{k}(:,steps_test),'k-','DisplayName','Full DMD Recon.')
        title(['Local HAVOK Forecast at t = ' num2str(t_test(1)) ': Component ' num2str(k)]);
        hold on
        plot(tHAVloc,xHAVloc,'r-','DisplayName','HAVOK Forecast')
        hold on
        plot([t_test(wSteps) t_test(wSteps)],ylim,'k--','DisplayName','Boundary Between Window and Future')
        hold off
        xlim([t_test(1) t_test(end)])
        ylim(yBounds)
    end
end
