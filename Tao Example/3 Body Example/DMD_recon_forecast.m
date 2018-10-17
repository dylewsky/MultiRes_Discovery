function xr = DMD_recon_forecast(t_test,wSteps,nComponents,recon_filter_sd,wHist,omegaHist,bHist,cHist,omClassHist,tStartHist,suppress_growth)
    % Given mwDMD results for nHist previous steps, forecast the state over
    % the duration of t_test (using a half-Gaussian interpolated filtering
    % scheme to weight historical contributions)
    %
    % t_test = time vector (first wSteps represent last observed DMD
    % window, the rest is the "future")
    %
    % nComponents = number of time scale classification components being
    % (separately) reconstructed
    % 
    % wSteps = width (in steps) of one window in moving-window DMD scheme
    % 
    % recon_filter_sd = standard deviation used in building Gaussian
    % filters for weighting of historical contributions
    %
    % [var]Hist are the values of DMD results for each parameter over the
    % last nHist windows
    %
    % suppressGrowth toggles the use of frequencies with real components >
    % 0 in the reconstruction

    nHist = size(omegaHist,1);
    nVars = size(omegaHist,2);
    if any([size(wHist,1), size(omegaHist,1), size(bHist,1), size(cHist,1), size(omClassHist,1), size(tStartHist,1)] - nHist ~= 0)
        disp('Input argument dimensional error')
        return;
    end
    if any([size(wHist,2), size(omegaHist,2), size(bHist,2), size(cHist,2), size(omClassHist,2)] - nVars ~= 0)
        disp('Input argument dimensional error')
        return;
    end
    
    t_test = t_test(:).'; %row vector
    dt = t_test(2) - t_test(1);
    stepSize = round((tStartHist(2) - tStartHist(1))/dt);
    t_test_pad = [(t_test(1) + (dt * (-wSteps:-1))) , t_test , (t_test(end) + (dt * (1:wSteps)))]; %pad w/ 2 extra windows worth of space
    
    if suppress_growth == 1
        omegaHist(real(omegaHist) > 0) = sqrt(-1) * imag(omegaHist(real(omegaHist) > 0));
    end
    
    
    % Build filters
    hFilter = exp(-((-stepSize + -wSteps+1 : (length(t_test)+wSteps)) - (3*wSteps/2)).^2/recon_filter_sd^2);
    % gaussian peaked at wSteps/2 (or 3*wSteps/2 if you count the
    % preliminary padding window)
    allFilters = zeros(nHist+1,length(t_test_pad));
    for ih = 1:nHist+1 % 1 extra element that will dominate at large t
        hShift = (ih-1)*stepSize; %extra element gets negative shift
        allFilters(ih,1 : end - hShift) = hFilter(hShift + 1 : end-stepSize);
    end
    pct_contributions = allFilters./(repmat(sum(allFilters),nHist+1,1));

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

    pct_contributions(:,end) = min(min(pct_contributions)) *ones(nHist,1);
    pct_contributions = fillmissing(pct_contributions.','pchip').';

    %normalize again
    pct_contributions = pct_contributions./repmat(sum(pct_contributions,1),nHist,1);

    allFilters = pct_contributions;
    
    % Create reconstruction
        
    xr = cell(nComponents,1); %stores separate time series reconstructions
    for k = 1:nComponents
        xr_comp = zeros(nVars,length(t_test_pad));
%         xn = zeros(1,length(t_test)); %only need to keep track of normalization for t_test
        xn = zeros(1,length(t_test_pad));
        xHist = zeros(nHist,nVars,length(t_test_pad));
        for ih = 1:nHist
            hShift = (ih-1)*stepSize;

            w = squeeze(wHist(ih,:,:));
            b = bHist(ih,:).';
            Omega = omegaHist(ih,:).';
            om_class = omClassHist(ih,:).';
            c = cHist(ih,:).';
            t_start = tStartHist(ih);

            xr_hWind = w(:, om_class == k)*diag(b(om_class == k))*exp(Omega(om_class == k)*(t_test_pad-t_start));
            if k == 1 % constant shift gets put in LF recon
                xr_hWind = xr_hWind + repmat(c,1,length(t_test_pad));
            end
            xHist(ih,:,:) = xr_hWind;
%             plot(t_test_pad,xr_hWind(1,:))
%             hold on
            
%             thisFilter = zeros(size(t_test_pad));
%             thisFilter(wSteps-hShift+1 : wSteps-hShift+length(t_test)) = allFilters(ih,:);
            thisFilter = allFilters(ih,:);

            xr_hWind = xr_hWind .* repmat(thisFilter,nVars,1);

%             xr_comp(:,wSteps-hShift+1 : end) = xr_comp(:,wSteps-hShift+1 : end) + xr_hWind(:,1:end-wSteps+hShift);
            xr_comp = xr_comp + xr_hWind;
%             xr_comp(:,wSteps+1 : end) = xr_comp(:,wSteps+1 : end).*repmat(thisFilter,nVars,1);
            %no need to apply filter to first wSteps, those will be cut
            %off anyway
            xn = xn + thisFilter;
        end
        xr_comp = real(xr_comp./repmat(xn,nVars,1)); %normalize
        xr_comp = xr_comp(:,wSteps+1 : end-wSteps); %cut off initial padding steps
        
        downsample_pad = 1:25:length(t_test_pad);
        
%         figure
%         for ih = 1:nHist
%             plot(t_test_pad(downsample_pad),allFilters(ih,downsample_pad),'Color',[ih/nHist, 0 0])
%             hold on
%         end
%         title('All historical filters')
%         plot(t_test_pad(downsample_pad),sum(allFilters(:,downsample_pad),1),'g')
%         
%         figure
%         for ih = 1:nHist
%             plot(t_test_pad(downsample_pad),(allFilters(ih,downsample_pad).').*squeeze(xHist(ih,1,downsample_pad)),'Color',[ih/nHist, 0 0])
%             hold on
%         end
%         title('All filtered DMD results')
%         
%         figure
%         for ih = 1:nHist
%             plot(t_test_pad(downsample_pad),squeeze(xHist(ih,1,downsample_pad)),'Color',[ih/nHist, 0 0])
%             hold on
%         end
%         plot(t_test_pad(downsample_pad),sum(allFilters(:,downsample_pad).*squeeze(xHist(:,1,downsample_pad)),1),'Color','g')
%         title('All historical DMD results')
        
        xr{k} = xr_comp;
    end
end

