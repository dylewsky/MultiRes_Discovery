clear; close all; clc

addpath('altmany-export_fig-9ac0917');
addpath(genpath(fullfile('..','Optimized DMD','optdmd-master')));
load('raw_data.mat');

r = size(x,1); %rank to fit w/ optdmd
imode = 1; %parameter for optdmd code
%  imode = 1, fit full data, slower
%  imode = 2, fit data projected onto first r POD modes
%      or columns of varargin{2} (should be at least r
%      columns in varargin{2})

windows = floor(2.^(9:0.5:14));
nLevels = length(windows);
nVars = size(x,1);
nSteps = 2^15; %truncation size of full data

stepSize = 2^6;

nSlide = floor((nSteps-windows)/stepSize);


%% execute optDMD
mrw_res = cell(nLevels,max(nSlide));

for lv = 1:nLevels
    lv_wind = 1:windows(lv);
    corner_sharpness = 64; %higher = sharper corners
    lv_kern = tanh(corner_sharpness*lv_wind/windows(1)) - tanh(corner_sharpness*(lv_wind-windows(1))/windows(1)) - 1;
    for k = 1:nSlide(lv)
        disp([lv,k]);
        sampleStart = stepSize*(k-1) + 1;
        sampleSteps = sampleStart : sampleStart + windows(lv) - 1;
        xSample = x(:,sampleSteps);
        tSample = TimeSpan(sampleSteps);
        
        mrw_res{lv,k}.x = xSample;
        xSample = xSample.*repmat(lv_kern,nVars,1); %round off corners
        mrw_res{lv,k}.t = tSample;
        [w, e, b] = optdmd(xSample,tSample,r,imode);
        mrw_res{lv,k}.w = w;
        mrw_res{lv,k}.Omega = e;
        mrw_res{lv,k}.b = b;
    end
end


%% Cluster Frequencies
close all;
if exist('mrw_res','var') == 0
    load('mrw_res.mat');
end

mrw_nz = zeros(nLevels, max(nSlide)); % tracks which cells have a real run associated with them
for lv = 1:nLevels
    for k = 1:nSlide(lv)
        mrw_nz(lv,k) = 1;
    end
end

nComponents = 2;
% nBins = 64;

gmmList = cell(nLevels,1);
kMeansList = cell(nLevels,1);
for lv = 1:nLevels
% for lv = 1:1
    all_om = [];
    for k = 1:nSlide(lv)
        if mrw_nz(lv,k) == 0
            continue
        end
        Omega = mrw_res{lv,k}.Omega;
        all_om = [all_om; Omega];
    end

    all_om_sq = sort(conj(all_om) .* all_om);
    all_om_sq = all_om_sq(1:floor(0.99*length(all_om_sq))); %remove top 1% as outliers
    
    [idx,clustCents,clustSumD,clustDists] = kmeans(all_om_sq,nComponents);
    
    [~,sortInd] = sort(clustCents);
    kMeansList{lv}.idx = idx;
    kMeansList{lv}.clustCents = clustCents;
    kMeansList{lv}.clustSumD = clustSumD;
    kMeansList{lv}.clustDists = clustDists;
    kMeansList{lv}.sortInd = sortInd;

    mean_clust_dist = 0;
    for k = 1:nSlide
        if mrw_nz(lv,k) == 0
            continue
        end
        omega = mrw_res{lv,k}.Omega;
        om_sq = omega.*conj(omega);
        
        om_sq_dist_compare = abs(repmat(om_sq,1,nComponents) - repmat(clustCents.',length(om_sq),1));
        [om_sq_dist,om_sq_clust] = min(om_sq_dist_compare,[],2);
        
        om_class = sortInd(om_sq_clust);
        om_mean_dist = norm(om_sq_dist);
        mean_clust_dist = mean_clust_dist + om_mean_dist;
        mrw_res{lv,k}.om_class = om_class;
        mrw_res{lv,k}.om_post = om_mean_dist;
    end
    mean_clust_dist = mean_clust_dist/nSlide(lv);
    kMeansList{lv}.mean_clust_dist = mean_clust_dist;
end

save('mrw_res.mat', 'mrw_res');    


%% Plot MultiRes Results
close all;
if exist('mrw_res','var') == 0
    load('mrw_res.mat');
end

export_result = 0;
logScale = 0;

% figure('units','pixels','Position',[0 0 1366 2*768])

% plotDims = [3 4]; %rows, columns of plot grid on screen at a given time
% plotDims = [1 4]; %rows, columns of plot grid on screen at a given time
colorList = {'b','r','g','k','y'};
    
x_PoT = x(:,1:nSteps);
t_PoT = TimeSpan(1:nSteps);
%res_list: [pn, level #, nSplit, sampleSteps/nSplit]

dupShift = 0.02; %multiplier to shift duplicate values so they can be visually distinguished

nBins = 64;

for lv = 1:nLevels
    figure('units','pixels','Position',[100 100 1200 400])
%     j = res_list(q,2);
%     pn = res_list(q,1);
%     nSplit = 2^(j-1);
    
    steps_per_window = windows(lv);
%     om_spec = zeros(nVars,nSteps);
    subplot(2,2,4);
    plot(t_PoT,real(x_PoT)) %plot ground truth
    xlabel('t')
    ylabel('Measurement Data')
    
    xMax = max(max(abs(x_PoT)));
    
    ylim(1.5*[-xMax, xMax]);
    hold on
    
    all_om = [];
    
    for k = 1:nSlide(lv)
        Omega = mrw_res{lv,k}.Omega;
        all_om = [all_om; Omega];
    end
    all_om_sq = sort(all_om .* conj(all_om));
    thresh_om_sq = all_om_sq(floor(0.99*length(all_om_sq)));
    all_om_sq = all_om_sq(1:floor(0.99*length(all_om_sq))); % cut off top 1%
    subplot(2,2,[1 3]);
    om_hist = histogram(all_om_sq,nBins);
    mesh_pad = 10;
    bin_mesh = om_hist.BinEdges(1)-mesh_pad:0.5:om_hist.BinEdges(end)+mesh_pad;
    xlabel('|\omega|^2');
    ylabel('Count');
    title('|\omega|^2 Spectrum & k-Means Centroids');
    hold on
    
    km_res = sort(kMeansList{lv}.clustCents);
    yl = ylim;
    for c = 1:length(km_res)
        plot([km_res(c), km_res(c)], yl,'Color',colorList{c},'LineWidth',2)
    end

    
    for k = 1:nSlide(lv)
        w = mrw_res{lv,k}.w;
%         e = mr_res{pn,j,k}.e;
        b = mrw_res{lv,k}.b;
        Omega = mrw_res{lv,k}.Omega;
        om_class = mrw_res{lv,k}.om_class;
        t = mrw_res{lv,k}.t;
        tShift = t-t(1); %compute each segment of xr starting at "t = 0"
        t_nudge = 5;
        
%         rankDefFlag = 0;
%         for bi = 1:length(b)
%             if b(bi) == 0
%                 w(:,bi) = zeros(size(w,1),1);
%                 rankDefFlag = 1;
%             end
%         end
        
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
        om_window_spec = repmat(om_sq, 1, steps_per_window);
%         om_spec(:,stepSize*(k-1) + 1 : stepSize*(k-1) + steps_per_window) = om_window_spec;

        subplot(2,2,2);
        
        for g = 1:nComponents
            om_window_spec_cat = om_window_spec(om_class == g,:);
            if isempty(om_window_spec_cat) == 1
                continue
            end
            if logScale == 1
                semilogy(t,om_window_spec_cat,'Color',colorList{g},'LineWidth',1);
            else
                plot(t,om_window_spec_cat,'Color',colorList{g},'LineWidth',1);
            end
            hold on
        end
        title('|\omega|^2 Spectra (Moving Window)');
        
%         xr_window = w*diag(b)*exp(Omega*t);
% 
%         
%         subplot(plotDims(1),plotDims(2),plotDims(2)-1:plotDims(2));
% %         if rankDefFlag == 0
% %             p_xr = plot(t,real(xr_window),'LineWidth',2);
% %         else
% %             p_xr = plot(t,real(xr_window),'-.','LineWidth',1);
% %         end
%         p_xr = plot(t,real(xr_window),'LineWidth',2);
%         title([num2str(steps_per_window) '-Step Window (Frequencies ~' num2str(1/(steps_per_window*(t(2)-t(1)))) ' Hz)']);
%         xlabel('t')
%         xlim([t_PoT(1) t_PoT(end)]);
%         ylabel('Re[x]')
%         hold on

    end
    
    if export_result == 1
        if lv == 1
            export_fig 'manyres_opt' '-pdf';
%             print(gcf, '-dpdf', 'manyres_opt.pdf'); 
        else
            export_fig 'manyres_opt' '-pdf' '-append';
%             print(gcf, '-dpdf', 'manyres_opt.pdf', '-append'); 
        end
        close(gcf);
    end
end
