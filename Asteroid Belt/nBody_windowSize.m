clear; close all; clc

% addpath('altmany-export_fig-9ac0917');
addpath('optdmd-master');
load('nBody_data.mat');

n = size(z,2)/4;
x = z(:,1:2*n).';
TimeSpan = t;


r = 8; %rank to fit w/ optdmd
imode = 1;
%  imode = 1, fit full data, slower
%  imode = 2, fit data projected onto first r POD modes
%      or columns of varargin{2} (should be at least r
%      columns in varargin{2})
nLevels = 5;
nComponents = 2;
nVars = size(x,1);
nSteps = 2^12;
downScale = 5; %different-sized data sets will be built from blocks of size nSteps/2^downScale

if mod(nSteps * 2^(-downScale) * 2^(-(nLevels-1)),2) ~= 0
    print('Error: nSteps not sufficiently divisible by 2')
    return;
end

primeList = primes(2^downScale);
primeList = [primeList(2:end) 2^downScale]; %remove superfluous 2 and add 2^downScale

primeList = primeList(end-4:end); %don't need to go to smallest sample sizes w/ such a long time series input

%% OptDMD

mr_res = cell(length(primeList),nLevels,2^(nLevels-1));
res_list = [];
for pn = 1:length(primeList)
    sampleSteps = nSteps * primeList(pn) / 2^(downScale);
    xSample = x(:,1:sampleSteps);
    tSample = TimeSpan(1:sampleSteps);
    
    nHold = 0;

    for n = (1:nLevels)+3
        nSplit = 2^(n-1);
        xL = reshape(xSample, nVars, sampleSteps/nSplit, nSplit);
        tL = reshape(tSample, 1, sampleSteps/nSplit, nSplit);
        
        res_list = [res_list; pn, n, nSplit, sampleSteps/nSplit];
        for k = 1:nSplit
            xT = xL(:,:,k);
            tT = tL(:,:,k);
            mr_res{pn,n,k}.x = xT;
            mr_res{pn,n,k}.t = tT;
            try
                [w, e, b] = optdmd(xT,tT,r,imode);
            catch ME
                disp(['Error on pn = ' num2str(pn) ', n = ' num2str(n) ', k = ' num2str(k)]);
                continue
            end
%             Omega = log(diag(e))/(TimeSpan(2) - TimeSpan(1));
            mr_res{pn,n,k}.w = w;
            mr_res{pn,n,k}.Omega = e;
            mr_res{pn,n,k}.b = b;
%             mr_res{pn,n,k}.Omega = Omega;
        end
    %     nHold = input('Subtract off how many modes?')
    end
end

% mr_res_stack = reshape(mr_res,length(primeList)*nLevels*2^(nLevels-1),1,1);
res_list = sortrows(res_list,4,'descend');


%% Cluster Frequencies
close all;
if exist('mr_res','var') == 0
    load('mr_res_2.mat');
    load('res_list_2.mat');
end

thresh_pct = 1;

for q = 1:size(res_list,1)
% for q = 16
    j = res_list(q,2);
    pn = res_list(q,1);
    nSplit = 2^(j-1);
    sampleSteps = nSteps * primeList(pn) / 2^(downScale);
    steps_per_window = sampleSteps/nSplit;
    om_spec = zeros(nVars,sampleSteps);
    nSlide = res_list(q,4);
    
    all_om = [];
    all_om_grouped = [];
    for k = 1:nSlide
        try
            Omega = mr_res{pn,n,k}.Omega;
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
            omega = mr_res{pn,n,k}.Omega;
        catch ME
            continue
        end
        om_sq = omega.*conj(omega);

        om_sq_dists  = (repmat(km_centroids.',r,1) - repmat(om_sq,1,nComponents)).^2;
        [~,om_class] = min(om_sq_dists,[],2);

        mr_res{pn,n,k}.om_class = om_class;
    end
    
end
save('mr_res_2.mat', 'mr_res');
save('res_list_2.mat', 'res_list');
save('gmm_list_2.mat', 'gmmList');

%% Plot MultiRes Results
close all;

export_result = 0;
logScale = 0;

if ~exist('res_list','var')
    load('res_list_2.mat');
    load('gmm_list_2.mat');
    load('mr_res_2.mat');
end

% figure('units','pixels','Position',[0 0 1366 2*768])

% plotDims = [3 4]; %rows, columns of plot grid on screen at a given time
plotDims = [1 4]; %rows, columns of plot grid on screen at a given time
colorList = {'b','r','g','k','y'};
    
x_PoT = x(:,1:nSteps);
t_PoT = TimeSpan(1:nSteps);
%res_list: [pn, level #, nSplit, sampleSteps/nSplit]

dupShift = 0.02; %multiplier to shift duplicate values so they can be visually distinguished

nBins = 64;

kmList = cell(size(res_list,2),1);

for q = 1:size(res_list,1)
% for q = [16]
    figure('units','pixels','Position',[100 100 1200 400])
    j = res_list(q,2);
    pn = res_list(q,1);
    nSplit = 2^(j-1);
    sampleSteps = nSteps * primeList(pn) / 2^(downScale);
    steps_per_window = sampleSteps/nSplit;
    om_spec = zeros(nVars,sampleSteps);
    omIm_spec = zeros(nVars,sampleSteps);
    b_spec = zeros(nVars,sampleSteps);
%     scrollsubplot(plotDims(1),plotDims(2),[plotDims(2)*q-1, plotDims(2)*q]);
    subplot(plotDims(1),plotDims(2),[plotDims(2)-1, plotDims(2)]);
    plot(t_PoT,real(x_PoT),'k-') %plot ground truth
    xMax = max(max(abs(x_PoT)));
    
    ylim(1.5*[-xMax, xMax]);
    hold on
    
    all_om = [];
    
    for k = 1:nSplit
        try
            Omega = mr_res{pn,j,k}.Omega;
        catch ME
            continue
        end
        all_om = [all_om; Omega];
    end
    all_om_sq = all_om .* conj(all_om);
    
    all_om_sq = sort(all_om_sq);
    all_om_sq_thresh = all_om_sq(floor(0.99*length(all_om_sq)));
    all_om_sq = all_om_sq(1:floor(0.99*length(all_om_sq)));
    
    subplot(plotDims(1),plotDims(2),plotDims(2)-2);
    om_hist = histogram(all_om_sq,nBins);
    mesh_pad = 10;
    bin_mesh = om_hist.BinEdges(1)-mesh_pad:0.5:om_hist.BinEdges(end)+mesh_pad;
    xlabel('|\omega|^2');
    ylabel('Count');
    hold on

    
    for k = 1:nSplit
        try
            w = mr_res{pn,j,k}.w;
        catch ME
            continue
        end
%         e = mr_res{pn,j,k}.e;
        b = mr_res{pn,j,k}.b;
        Omega = mr_res{pn,j,k}.Omega;
        
%         if isempty(gmmList{q}) == 0
%             om_class = mr_res{pn,j,k}.om_class;
%         end
        try
            om_class = mr_res{pn,j,k}.om_class;
        catch ME
            continue
        end
        t = mr_res{pn,j,k}.t;
        tShift = t-t(1); %compute each segment of xr starting at "t = 0"
        t_nudge = 5;
        
        rankDefFlag = 0;
        for bi = 1:length(b)
            if b(bi) == 0
                w(:,bi) = zeros(size(w,1),1);
                rankDefFlag = 1;
            end
        end
        
        % Plot |omega|^2 spectrum
        om_sq = conj(Omega).*Omega;
        for oi = 1:length(om_sq)
            for oj = oi:length(om_sq)
                if (abs(om_sq(oi)-om_sq(oj))/((om_sq(oi)+om_sq(oj))/2)) <= dupShift && (oi ~= oj)
                    om_sq(oi) = om_sq(oi)*(1-dupShift);
                    om_sq(oj) = om_sq(oj)*(1+dupShift);
                end
            end
        end
        om_window_spec = repmat(om_sq, 1, steps_per_window);
        om_spec(:,(k-1)*steps_per_window+1:k*steps_per_window) = om_window_spec;
        
        subplot(plotDims(1),plotDims(2),plotDims(2)-3);
        if isempty(gmmList{q}) == 0 %plot clustered frequencies if GMM exists
            for g = 1:nComponents
                om_window_spec_cat = om_window_spec(om_class == g,:);
                if isempty(om_window_spec_cat) == 1
                    continue
                end
                if logScale == 1
                    semilogy(t,om_window_spec_cat,'Color',colorList{g},'LineWidth',3);
                else
                    plot(t,om_window_spec_cat,'Color',colorList{g},'LineWidth',3);
                end
                hold on
            end
        else %if there weren't enough data points to train GMM, just plot spectrum
            if logScale == 1
                semilogy(t,om_window_spec,'LineWidth',3);
            else
                plot(t,om_window_spec,'LineWidth',3);
            end
        end
        ylim([0 all_om_sq_thresh]);
%         p_om_sq = semilogy(t,om_window_spec,'LineWidth',3);
%         title(['Frequency Spectrum for ' num2str(steps_per_window) '-Step Window']);
        xlabel('t')
        xlim([t_PoT(1) t_PoT(end)]);
        ylabel('| \omega |^2')
        hold on
        
%         % Plot |Im[omega]| spectrum
%         omIm = abs(imag(Omega));
%         for oi = 1:length(omIm)
%             for oj = oi:length(omIm)
%                 if (abs(omIm(oi)-omIm(oj))/((omIm(oi)+omIm(oj))/2)) <= dupShift && (oi ~= oj)
%                     omIm(oi) = omIm(oi)*(1-dupShift);
%                     omIm(oj) = omIm(oj)*(1+dupShift);
%                 end
%             end
%         end
%         omIm_window_spec = repmat(omIm, 1, steps_per_window);
%         omIm_spec(:,(k-1)*steps_per_window+1:k*steps_per_window) = omIm_window_spec;
%         
%         subplot(plotDims(1),plotDims(2),plotDims(2)-2);
%         p_om_im = semilogy(t,omIm_window_spec,'LineWidth',3);
% %         title(['Frequency Spectrum for ' num2str(steps_per_window) '-Step Window']);
%         xlabel('t')
%         xlim([t_PoT(1) t_PoT(end)]);
%         ylabel('|Im[omega]|')
%         hold on

%         if j ~=nLevels
%             set(gca,'XTick',[])
%         end
        
      
        
%         if j ~=nLevels
%             set(gca,'XTick',[])
%         end
        
%         xr_window = zeros(nVars,steps_per_window);
%         for m = 1:nVars
%             xr_window = xr_window + Phi(:,m) * exp(w(m)*tShift) * b(m);
%         end
        
%         xr_window = w*diag(b)*exp(Omega*tShift);
        xr_window = w*diag(b)*exp(Omega*t);

        
%         scrollsubplot(plotDims(1),plotDims(2),plotDims(2)*q-1:plotDims(2)*q);
        subplot(plotDims(1),plotDims(2),plotDims(2)-1:plotDims(2));
        if rankDefFlag == 0
            p_xr = plot(t,real(xr_window),'LineWidth',2);
        else
            p_xr = plot(t,real(xr_window),'-.','LineWidth',1);
        end
        title([num2str(steps_per_window) '-Step Window (Frequencies ~' num2str(1/(steps_per_window*(t(2)-t(1)))) ' Hz)']);
        xlabel('t')
        xlim([t_PoT(1) t_PoT(end)]);
        ylabel('Re[x]')
        hold on

%         if j ~=nLevels
%             set(gca,'XTick',[])
%         end

% %         b = b ./ xr_window(:,1); %normalize so weights are as though the window starts at t=0
%         b_sq = conj(b).*b;
%         b_window_spec = repmat(b_sq, 1, steps_per_window);
% %         scrollsubplot(plotDims(1),plotDims(2),plotDims(2)*q-2);
%         subplot(plotDims(1),plotDims(2),plotDims(2)-2);
%         plot(t,b_window_spec,'LineWidth',2);
% %         title(['Weights for ' num2str(steps_per_window) '-Step Window']);
%         xlabel('t')
%         xlim([t_PoT(1) t_PoT(end)]);
%         ylim auto;
%         ylabel('|b|^2')
%         hold on
    end
    for k = 1:nSplit %plot dotted lines between time splits
        t = mr_res{pn,j,k}.t;
%         scrollsubplot(plotDims(1),plotDims(2),plotDims(2)*q-3);
        subplot(plotDims(1),plotDims(2),plotDims(2)-3);
    	plot([t(end) t(end)],get(gca, 'YLim'),'k:')
        hold on 
% %         scrollsubplot(plotDims(1),plotDims(2),plotDims(2)*q-2);
%         subplot(plotDims(1),plotDims(2),plotDims(2)-2);
%         plot([t(end) t(end)],get(gca, 'YLim'),'k:')
%         hold on 
%         scrollsubplot(plotDims(1),plotDims(2),plotDims(2)*q-1:plotDims(2)*q);
        subplot(plotDims(1),plotDims(2),plotDims(2)-1:plotDims(2));
        plot([t(end) t(end)],get(gca, 'YLim'),'k:')
        hold on 
    end
    if export_result == 1
        if q == 1
            export_fig 'manyres_opt_2' '-pdf';
%             print(gcf, '-dpdf', 'manyres_opt.pdf'); 
        else
            export_fig 'manyres_opt_2' '-pdf' '-append';
%             print(gcf, '-dpdf', 'manyres_opt.pdf', '-append'); 
        end
        close(gcf);
    end
end
