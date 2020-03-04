%%% DMD SCALE SEPARATION FOR A THE ORBITAL SYSTEM OF JUPITER, SATURN, AND THE SUN
%%%
%%% Before running, make sure to 1) download and set up the optdmd package
%%% available here: https://github.com/duqbo/optdmd and 2) run the
%%% simulation Three_Body_Orbit_Generate_Data.m to generate the data
%%%
%%% This code reproduces the results of Sec. V of "Dynamic mode decomposition
%%% for multiscale nonlinear physics," by Dylewsky, Tao, and Kutz
%%% (DOI: 10.1103/PhysRevE.99.063311)


clear variables; close all; clc

colorList = {'k','r','b','g','c','m'};

imode = 1; %parameter for optdmd code
%  imode = 1, fit full data, slower
%  imode = 2, fit data projected onto first r POD modes
%      or columns of varargin{2} (should be at least r
%      columns in varargin{2})
global_SVD = 1; %use global SVD modes for each DMD rather than individual SVD on each window

%global settings for recursions with initialize_artificially == 1
use_last_freq = 0; % initialize optDMD with eigenvalues from previous window
use_median_freqs = 1; % initialize optDMD with median eigenvalues of previous run
    
n_recursions = 4;

rec_x = {'Three_Body_Data_Cartesian.mat',...
    'mwDMD_sep_recon_recursion_01.mat',...
    'mwDMD_sep_recon_recursion_02.mat',...
    'mwDMD_sep_recon_recursion_03.mat'};
rec_composite_components = {1, 1, 1, []}; %scale components which get subdivided in subsequent recursions
rec_wSteps = [6000, 460, 7200, 19200]; % window sizes for each recursion iteration
rec_stepSize = [400, 20, 80, 160];
rec_nComponents = [3, 2, 2, 2]; % Number of components to seek in each recursion iteration
rec_constrain_LHP = [1, 1, 1, 1]; % pass constraint to optDMD to force eigenvalues into left half of complex plane
rec_r = 2*rec_nComponents;
rec_initialize_artificially = [1 0 0 0];

save('recursive_params.mat','n_recursions','rec_x','rec_composite_components','rec_wSteps',...
    'rec_stepSize','rec_nComponents','rec_r','rec_initialize_artificially','rec_constrain_LHP');

%% Run Recursive Sliding-Window DMD
for nr = 1:n_recursions
    load(rec_x{nr}); % Results are saved between each recursion
    if nr == 1
        x = pos;
    else
        x = xr_sep{rec_composite_components{nr-1}}; %use redundant component of previous recursion
    end
    TimeSpan = tspan;

    wSteps = rec_wSteps(nr);
    stepSize = rec_stepSize(nr);
    nComponents = rec_nComponents(nr);
    constrain_LHP = rec_constrain_LHP(nr);
    r = rec_r(nr);
    initialize_artificially = rec_initialize_artificially(nr);

    if global_SVD == 1
        [u,~,v] = svd(x,'econ');
        u = u(:,1:r); %just first r modes
        v = v(:,1:r);
    end
    
    if initialize_artificially == 1
        if use_median_freqs == 1
            load(['km_centroids_recursion_' num2str(nr,'%02.f') '.mat']);
            freq_meds = repelem(sqrt(km_centroids),r/nComponents);
        end
    end
    
    nSplit = floor(length(TimeSpan)/wSteps); %number of windows if they were non-overlapping
    nSteps = wSteps * nSplit;
    nVars = size(x,1);
    thresh_pct = 1;

    nSlide = floor((nSteps-wSteps)/stepSize);

    save(['mwDMD_params_recursion_' num2str(nr,'%02.f') '.mat'],'r','nComponents','initialize_artificially','use_last_freq','use_median_freqs','wSteps','nSplit','nSteps','nVars','thresh_pct','stepSize','nSlide');

    %% execute optDMD

    corner_sharpness = 16; %higher = sharper corners
    lv_kern = tanh(corner_sharpness*(1:wSteps)/wSteps) - tanh(corner_sharpness*((1:wSteps)-wSteps)/wSteps) - 1;

    mr_res = cell(nSlide,1);
    t_starts = zeros(nSlide,1);
    for k = 1:nSlide
        sampleStart = stepSize*(k-1) + 1;
        sampleSteps = sampleStart : sampleStart + wSteps - 1;
        xSample = x(:,sampleSteps);
        tSample = TimeSpan(sampleSteps);

        mr_res{k}.x = xSample;
        mr_res{k}.t = tSample;

        c = mean(xSample,2); %subtract off mean before rounding corners
        xSample = xSample - repmat(c,1,size(xSample,2));

        xSample = xSample.*repmat(lv_kern,nVars,1); %round off corners

        t_start = tSample(1);
        tSample = tSample - t_start;
        if global_SVD == 0
            [u,~,~] = svd(xSample,'econ');
            u = u(:,1:r);
        end
        
        % Constrain eigenvalues (if applicable)
        lbc = [-Inf*ones(r,1); -Inf*ones(r,1)];
        if constrain_LHP == 1
            ubc = [zeros(r,1); Inf*ones(r,1)];
        else
            ubc = [Inf*ones(r,1); Inf*ones(r,1)];
        end
        copts = varpro_lsqlinopts('lbc',lbc,'ubc',ubc);
        
        if (exist('e_init','var')) && (initialize_artificially == 1)
            try
                [w, e, b] = optdmd(xSample,tSample,r,imode,[],e_init,u,copts);
            catch ME
                run('setup.m');
                [w, e, b] = optdmd(xSample,tSample,r,imode,[],e_init,u,copts);
            end

        else
            try
                [w, e, b] = optdmd(xSample,tSample,r,imode,[],[],u,copts);
            catch ME
                run('setup.m');
                [w, e, b] = optdmd(xSample,tSample,r,imode,[],[],u,copts);
            end
        end
        if initialize_artificially == 1
            if use_last_freq == 1
                e_init = e;
            end
            if use_median_freqs == 1
                [eSq, eInd] = sort(e.*conj(e)); %match order to that of freq_meds
                freq_angs = angle(e(eInd));
                e_init = exp(sqrt(-1)*freq_angs) .* freq_meds;
            end
        end
        mr_res{k}.w = w;
        mr_res{k}.Omega = e;
        mr_res{k}.b = b;
        mr_res{k}.c = c;
        mr_res{k}.t_start = t_start;
        t_starts(k) = t_start;
    end


    %% Cluster Frequencies
    close all;
    if exist('mr_res','var') == 0
        try
            load(['mwDMD_mr_res_recursion_' num2str(nr,'%02.f') '_i2.mat']);
        catch ME
            load(['mwDMD_mr_recursion_' num2str(nr,'%02.f') '_res.mat']);
        end
    end


    all_om = [];
    all_om_grouped = [];
    for k = 1:nSlide
        try
            Omega = mr_res{k}.Omega;
        catch ME
            disp(['Error on k = ' num2str(k)]);
            continue
        end
        all_om = [all_om; Omega];
        [~,imSortInd] = sort(imag(Omega));
        all_om_grouped = [all_om_grouped; Omega(imSortInd).'];
    end

    all_om_sq = conj(all_om) .* all_om;
    all_om_sq = sort(all_om_sq);
    all_om_sq_thresh = all_om_sq(floor(thresh_pct*length(all_om_sq)));
    all_om_sq = all_om_sq(1:floor(thresh_pct*length(all_om_sq)));


    [~, km_centroids] = kmeans(all_om_sq,nComponents,'Distance','cityblock','Replicates',5);

    [km_centroids,sortInd] = sort(km_centroids);

    for k = 1:nSlide
        try
            omega = mr_res{k}.Omega;
        catch ME
            continue
        end
        om_sq = omega.*conj(omega);

        om_sq_dists  = (repmat(km_centroids.',r,1) - repmat(om_sq,1,nComponents)).^2;
        [~,om_class] = min(om_sq_dists,[],2);

        mr_res{k}.om_class = om_class;
    end

    %% Save mr_res

    if initialize_artificially == 0
        save(['mwDMD_mr_res_recursion_' num2str(nr,'%02.f') '.mat'], 'mr_res', 't_starts', '-v7.3');
    else
        save(['mwDMD_mr_res_recursion_' num2str(nr,'%02.f') '_i2.mat'], 'mr_res', 't_starts', '-v7.3');
    end


    %% Plot MultiRes Results
    close all;
    if exist('mr_res','var') == 0
        try
            load(['mwDMD_mr_res_recursion_' num2str(nr,'%02.f') '_i2.mat']);
        catch ME
            load(['mwDMD_mr_res_recursion_' num2str(nr,'%02.f') '.mat']);
        end
    end

    logScale = 0;

    x_PoT = x(:,1:nSteps);
    t_PoT = TimeSpan(1:nSteps);

    dupShift = 0; %multiplier to shift duplicate values so they can be visually distinguished

    nBins = 64;

    figure('units','pixels','Position',[100 100 1200 400])

    all_om = [];


    for k = 1:nSlide
        Omega = mr_res{k}.Omega;
        all_om = [all_om; Omega];
    end
    all_om_sq = sort(all_om .* conj(all_om));
    thresh_om_sq = all_om_sq(floor(thresh_pct*length(all_om_sq)));
    all_om_sq = all_om_sq(1:floor(thresh_pct*length(all_om_sq)));
    subplot(2,2,[1 3]);
    om_hist = histogram(all_om_sq,nBins);
    mesh_pad = 10;
    bin_mesh = om_hist.BinEdges(1)-mesh_pad:0.5:om_hist.BinEdges(end)+mesh_pad;
    xlabel('|\omega|^2');
    ylabel('Count');
    title('|\omega|^2 Spectrum & k-Means Centroids');
    hold on

    yl = ylim;
    for c = 1:nComponents
        plot([km_centroids(c), km_centroids(c)], yl,'Color',colorList{c},'LineWidth',2)
    end
    xr = zeros(size(x(:,1:nSteps))); %reconstructed time series
    xn = zeros(nSteps,1); %count # of windows contributing to each step
    for k = 1:nSlide
        w = mr_res{k}.w;
        b = mr_res{k}.b;
        Omega = mr_res{k}.Omega;
        om_class = mr_res{k}.om_class;
        t = mr_res{k}.t;
        c = mr_res{k}.c;
        t_start = mr_res{k}.t_start;
        tShift = t-t(1); %compute each segment of xr starting at "t = 0"
        xr_window = w*diag(b)*exp(Omega*(t-t_start)) + c;
        xr(:,(k-1)*stepSize+1:(k-1)*stepSize+wSteps) = xr(:,(k-1)*stepSize+1:(k-1)*stepSize+wSteps) + xr_window;
        xn((k-1)*stepSize+1:(k-1)*stepSize+wSteps) = xn((k-1)*stepSize+1:(k-1)*stepSize+wSteps) + 1;

        % Plot |omega|^2 spectrum
        om_sq = conj(Omega).*Omega;

        om_class = om_class(om_sq <= thresh_om_sq);
        om_sq = om_sq(om_sq <= thresh_om_sq); %threshold to remove outliers
        for oi = 1:length(om_sq)
            for oj = oi:length(om_sq)
                if (abs(om_sq(oi)-om_sq(oj))/((om_sq(oi)+om_sq(oj))/2)) <= dupShift && (oi ~= oj)
                    om_sq(oi) = om_sq(oi)*(1-dupShift);
                    om_sq(oj) = om_sq(oj)*(1+dupShift);
                end
            end
        end

        subplot(2,2,2);

        for g = 1:nComponents
            om_window_spec_cat = om_sq(om_class == g,:);
            if isempty(om_window_spec_cat) == 1
                continue
            end
            if logScale == 1
                for f = 1:length(om_window_spec_cat)
                    semilogy(mean(t),om_window_spec_cat(f),'.','Color',colorList{g},'LineWidth',1,'MarkerSize',10);
                end
            else
                for f = 1:length(om_window_spec_cat)
                    plot(mean(t),om_window_spec_cat,'.','Color',colorList{g},'LineWidth',1,'MarkerSize',10);
                end
            end
            hold on
        end
        title('|\omega|^2 Spectra (Moving Window)');
    end

    subplot(2,2,4);
    plot(t_PoT,real(x_PoT),'k-','LineWidth',1.5) %plot ground truth
    xlabel('t')
    ylabel('Measurement Data')
    xMax = max(max(abs(x_PoT)));
    ylim(1.5*[-xMax, xMax]);
    hold on

    xr = xr./repmat(xn.',nVars,1); %weight xr so all steps are on equal footing
    plot(t_PoT,real(xr),'b-','LineWidth',1.5) %plot averaged reconstruction
    title('Input (Black); DMD Recon. (Blue)');


    %% Link Consecutive Modes
    if exist('mr_res','var') == 0
        try
            load(['mwDMD_mr_res_recursion_' num2str(nr,'%02.f') '_i2.mat']);
        catch ME
            load(['mwDMD_mr_res_recursion_' num2str(nr,'%02.f') '.mat']);
        end
    end

    allModes = zeros(nSlide,nVars,r);
    allFreqs = zeros(nSlide,r);
    catList = flipud(perms(1:r));
    modeCats = zeros(nSlide,r);
    modeCats(1,:) = 1:r; %label using order of 1st modes
    allModes(1,:,:) = mr_res{1}.w;
    allFreqs(1,:) = mr_res{1}.Omega;
    [~,catsAsc] = sort(mr_res{1}.Omega.*conj(mr_res{1}.Omega)); %category labels in ascending frequency order

    distThresh = 0; %if more than one mode is within this distance then match using derivative
    permDistsHist = zeros(size(catList,1),nSlide-1);
    for k = 2:nSlide
        w_old = squeeze(allModes(k-1,:,:));
        w = mr_res{k}.w;
        permDists = zeros(size(catList,1),1);
        for np = 1:size(catList,1)
            permDists(np) = norm(w(:,catList(np,:))-w_old);
        end
        permDistsHist(:,k-1) = sort(permDists);
    %     minDists(k-1) = min(permDists);
        if nnz(permDists <= distThresh) >= 1
            if k == 2 %can't compute derivative 1 step back, so skip it
                [~,minCat] = min(permDists);
                continue;
            end
            w_older = squeeze(allModes(k-2,:,:));
            [~,minInds] = sort(permDists); %minInds(1) and minInds(2) point to 2 lowest distances

            wPrime = w - w_old;
            wPrime_old = w_old - w_older;

            for j = 1:size(wPrime,2)
                wPrime(:,j) = wPrime(:,j)/norm(wPrime(:,j));
                wPrime_old(:,j) = wPrime_old(:,j)/norm(wPrime_old(:,j));
            end

            derivDists = zeros(2,1);
            for np = 1:length(derivDists)
                derivDists(np) = norm(wPrime(:,catList(minInds(np),:))-wPrime_old);
            end

            [~,minDerInd] = min(derivDists);
            derivDistsSort = sort(derivDists);
            permDistsSort= sort(permDists);
            disp(['Derivative-Matching on k = ' num2str(k)])
            if (derivDistsSort(2)-derivDistsSort(1)) > (permDistsSort(2) - permDistsSort(1))
                minCat = minInds(minDerInd);
                [~,minCatNaive] = min(permDists);
                if minCat ~= minCatNaive
                    disp('Naive result was wrong!')
                end
            else
                [~,minCat] = min(permDists);
                disp('Naive result trumps derivative match')
            end

        else
            [~,minCat] = min(permDists);
        end
        modeCats(k,:) = catList(minCat,:);
        allModes(k,:,:) = w(:,catList(minCat,:));
        allFreqs(k,:) = mr_res{k}.Omega(catList(minCat,:));
    end
    meanFreqsSq = sum(allFreqs.*conj(allFreqs),1)/nSlide;


    if initialize_artificially == 0
        save(['mwDMD_allModes_recursion_' num2str(nr,'%02.f') '.mat'],'allModes','allFreqs');
    else
        save(['mwDMD_allModes_recursion_' num2str(nr,'%02.f') '_i2.mat'],'allModes','allFreqs');    
    end

    %% Second derivative of mode coordinates
    sd_norms = zeros(r,size(allModes,1)-2);
    for j = 1:r
        tsPP = allModes(1:end-2,:,j) + allModes(3:end,:,j) - 2*allModes(2:end-1,:,j);
        sd_norms(j,:) = vecnorm(tsPP.');
    end
    sd_norms = vecnorm(sd_norms).';

    sd_thresh = 0.1; %threshold above which to try different permutations to reduce 2nd deriv
                     %can be set to 0 to just check every step
    ppOverrideRatio = 0.7; %apply new permutation if 2nd deriv is reduced by at least this factor

    for k = 2:nSlide-1
        thisData = allModes(k-1:k+1,:,:);
        tdPP = squeeze(thisData(1,:,:) + thisData(3,:,:) - 2*thisData(2,:,:));
        if norm(tdPP) < sd_thresh
            continue;
        end
        permDists = zeros(size(catList,1),1);
        permData = zeros(size(thisData));
        for np = 1:size(catList,1)
            %hold k-1 configuration constant, permute k and k+1 registers
            permData(1,:,:) = thisData(1,:,:);
            permData(2:3,:,:) = thisData(2:3,:,catList(np,:));
            permPP = squeeze(permData(1,:,:) + permData(3,:,:) - 2*permData(2,:,:));
            permDists(np) = norm(permPP);
            if permDists(np) <= ppOverrideRatio * norm(tdPP)
                newCat = catList(np,:);
                allModes(k:end,:,:) = allModes(k:end,:,newCat);
                modeCats(k,:) = modeCats(k,newCat);
                disp(['2nd Deriv. Perm. Override: k = ' num2str(k)]);
            end
        end
    end

    if initialize_artificially == 0
        save(['mwDMD_allModes_recursion_' num2str(nr,'%02.f') '.mat'],'allModes','allFreqs');
    else
        save(['mwDMD_allModes_recursion_' num2str(nr,'%02.f') '_i2.mat'],'allModes','allFreqs');    
    end

    %% Split Regime Reconstruction

    suppress_growth = 1; %kill positive real components of frequencies

    xr_sep = cell(nComponents,1);
    for j = 1:nComponents
        xr_sep{j} = zeros(size(x(:,1:nSteps)));
    end
    
    dt = mr_res{1}.t(2)-mr_res{1}.t(1);
    all_b = zeros(r,nSlide);
    xn = zeros(nSteps,1); %track total contribution from all windows to each time step
    % Convolve each windowed reconstruction with a gaussian
    % Note that this will leave the beginning & end of time series prone to a
    % lot of error. Best to cut off wSteps from each end before further
    % analysis
    recon_filter_sd = wSteps/8; %std dev of gaussian filter
    recon_filter = exp(-((1 : wSteps) - (wSteps+1)/2).^2/recon_filter_sd^2);

    for k = 1:nSlide
        w = mr_res{k}.w;
        b = mr_res{k}.b;
        all_b(:,k) = b;
        Omega = mr_res{k}.Omega;
        if suppress_growth == 1
            Omega(real(Omega) > 0) = sqrt(-1) * imag(Omega(real(Omega) > 0));
        end
        om_class = mr_res{k}.om_class;
        t = mr_res{k}.t;
        c = mr_res{k}.c;
        t_start = mr_res{k}.t_start;
        tShift = t-t(1); %compute each segment of xr starting at "t = 0"
    %     t_nudge = 5;
        xr_sep_window = cell(nComponents,1);
        for j = 1:nComponents
            xr_sep_window{j} = w(:, om_class == j)*diag(b(om_class == j))*exp(Omega(om_class == j)*(t-t_start));
            if j == 1 % constant shift gets put in LF recon
                xr_sep_window{j} = xr_sep_window{j} + c;
            end
            xr_sep_window{j} = xr_sep_window{j}.*repmat(recon_filter,nVars,1);
            xr_sep{j}(:,(k-1)*stepSize+1:(k-1)*stepSize+wSteps) = xr_sep{j}(:,(k-1)*stepSize+1:(k-1)*stepSize+wSteps) + xr_sep_window{j};
        end

        xn((k-1)*stepSize+1:(k-1)*stepSize+wSteps) = xn((k-1)*stepSize+1:(k-1)*stepSize+wSteps) + recon_filter.';
%         window_contributions((k-1)*stepSize+1:(k-1)*stepSize+wSteps,k) = recon_filter.';
    end
    for j = 1:nComponents
        xr_sep{j} = xr_sep{j}./repmat(xn.',nVars,1);
    end

    tspan = TimeSpan(1:nSteps);

    % figure
    % plot(all_b.')
    % title('b values')

    nanCols = isnan(vertcat(xr_sep{:}));
    nanCols = sum(nanCols,1) ~= 0;


    for j = 1:nComponents
        xr_sep{j} = real(xr_sep{j}(:,~nanCols));
    end
    tspan = tspan(~nanCols);
    
    if nr == 1 % On first iteration, downsample output to save on memory
        ds_factor = 100;

        xr_sep_ds = cell(size(xr_sep));
        for j = 1:length(xr_sep)
            xr_sep_ds{j} = xr_sep{j}(:,1:ds_factor:end);
        end
        tspan = tspan(1:ds_factor:end);
        xr_sep = xr_sep_ds;
    end
    
    save(['recon_gaussian_filter_recursion_' num2str(nr,'%02.f') '.mat'],'recon_filter','recon_filter_sd');
    save(['mwDMD_sep_recon_recursion_' num2str(nr,'%02.f') '.mat'],'xr_sep','tspan','nComponents','suppress_growth');

    %% Plot Separated Reconstruction
    figure
    subplot(nComponents+1,1,1)
    plot(tspan,x(:,1:length(tspan)))
    title('Input Data')
    for j = 1:nComponents
        subplot(nComponents+1,1,j+1)
        plot(tspan,xr_sep{j})
        title(['Component ' num2str(j) ' Reconstruction'])
    end

    % Some optional diagnostic plots:
    
    % %% SVD On Reconstructions
    % U_sep = cell(nComponents,1);
    % S_sep = U_sep;
    % V_sep = U_sep;
    % r_sep = [4 2 2]; %truncation ranks for each component
    % figure
    % for j = 1:nComponents
    %     [U_sep{j},S_sep{j},V_sep{j}] = svd(xr_sep{j},'econ');
    %     subplot(1,nComponents,j)
    %     plot(cumsum(diag(S_sep{j}))/sum(diag(S_sep{j})),'*-','LineWidth',2)
    %     hold on
    %     plot([r_sep(j) r_sep(j)],ylim,'k:')
    %     title(['SVD Spectrum: Component #' num2str(j)])
    %     U_sep{j} = U_sep{j}(:,1:r_sep(j));
    %     S_sep{j} = S_sep{j}(1:r_sep(j),1:r_sep(j));
    %     V_sep{j} = V_sep{j}(:,1:r_sep(j));
    % end
    % 
    % figure
    % subplot(nComponents+1,1,1)
    % plot(TimeSpan(1:nSteps),x(:,1:nSteps))
    % title('Input Data')
    % for j = 1:nComponents
    %     subplot(nComponents+1,1,j+1)
    %     plot(tspan,V_sep{j})
    %     title(['Component ' num2str(j) ' Reconstruction: Top ' num2str(r_sep(j)) ' Modes'])
    % end
    % 
    % %% Compare Slow Reconstruction to Orbital Parameters
    % load('Three_Body_Data_Orbital.mat');
    % plotDomain = [tspan(wSteps) tcoarse(end)-tspan(wSteps)];
    % figure
    % subplot(3,1,1)
    % plot(tspan,V_sep{1})
    % title(['Slow-Regime Reconstruction: Top ' num2str(r_sep(1)) ' SVD Modes'])
    % xlim(plotDomain)
    % subplot(3,1,2)
    % plot(tcoarse,OOmega - mean(OOmega,2), ...
    %     tcoarse,aa - mean(aa,2), ...
    %     tcoarse,ee - mean(ee,2), ...
    %     tcoarse,ii - mean(ii,2), ...
    %     tcoarse,oomega - mean(oomega,2))
    % title('Slow Orbital Parameters (Zero-Average)')
    % % legend
    % xlim(plotDomain)
    % subplot(3,1,3)
    % [U_orb,S_orb,V_orb] = svd([OOmega;aa;ee;ii;oomega],'econ');
    % plot(tcoarse,V_orb(:,1:r_sep(1)))
    % xlim(plotDomain)
    % title(['Slow Orbital Parameters: Top ' num2str(r_sep(1)) ' SVD Modes'])
    % % figure
    % % plot(cumsum(diag(S_orb))/sum(diag(S_orb)),'*-','LineWidth',2)
    % % title('Orbital Parameters SVD Spectrum')
    % clear('plotDomain')
    % 
    % %% Construct A Matrices
    % if global_SVD ~= 1
    %     disp('Warning: Averaging A matrices in different bases')
    % end
    % A_avg = zeros(r);
    % A_proj_avg = zeros(r);
    % A_thresh = 100;
    % for k = 1:nSlide
    %     x_window = mr_res{k}.x;
    %     w = mr_res{k}.w;
    %     Omega = mr_res{k}.Omega;
    %     A = w*diag(Omega)*pinv(w);
    %     A_proj =  (u'*w)*diag(Omega)*(pinv(w)*u);
    %     
    %     A(abs(A) <= A_thresh) = 0; %threshold small values
    %     A_proj(abs(A_proj) <= A_thresh) = 0;
    %     
    %     A_avg = A_avg + A;
    %     A_proj_avg = A_proj_avg + A_proj;
    % end
    % A_avg = real(A_avg)/nSlide; %strips small imaginary residuals
    % A_proj_avg = real(A_proj_avg)/nSlide;
    % 
    % figure
    % plot(TimeSpan(1:nSteps),v(1:nSteps,:))
    % legend({'Mode 1','Mode 2','Mode 3','Mode 4'})
    % 
    % %% Time series of mode projections
    % allB = zeros(nSlide,r);
    % colFreqs = zeros(r,1);
    % for k = 1:nSlide
    %     thisCat = modeCats(k,:);
    %     b = mr_res{k}.b;
    % %     b = b(thisCat);
    %     allB(k,:) = b.';
    %     colFreqs = colFreqs + abs(mr_res{k}.Omega);
    % end
    % colFreqs = colFreqs/nSlide;
    % wMids = (TimeSpan(2)-TimeSpan(1))*(wSteps/2 + stepSize*(1:nSlide));
    % figure
    % subplot(2,1,1)
    % plot(wMids,allB(:,[1 3]))
    % legend({'LF','HF'});
    % title('OptDMD Weight Coefficients Over Time');
    % xlim([wMids(1) wMids(end)]);
    % subplot(2,1,2)
    % plot(TimeSpan(1:nSteps),x(:,1:nSteps))
    % xlim([wMids(1) wMids(end)]);

end