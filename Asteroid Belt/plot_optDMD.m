function plot_optDMD(mr_res,km_centroids,x,TimeSpan,r,wSteps,nSplit,stepSize,nComponents,thresh_pct)
    nSteps = wSteps * nSplit;
    nVars = size(x,1);
    nSlide = floor((nSteps-wSteps)/stepSize);
    colorList = 0.9*rand(nComponents,3);

    logScale = 0;
    
    x_PoT = x(:,1:nSteps);
    t_PoT = TimeSpan(1:nSteps);
    %res_list: [pn, level #, nSplit, sampleSteps/nSplit]

    dupShift = 0; %multiplier to shift duplicate values so they can be visually distinguished

    nBins = 64;

    figure('units','pixels','Position',[100 100 1200 400])
    %     j = res_list(q,2);
    %     pn = res_list(q,1);
    %     nSplit = 2^(j-1);

    %     om_spec = zeros(nVars,nSteps);
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
%     mesh_pad = 10;
%     bin_mesh = om_hist.BinEdges(1)-mesh_pad:0.5:om_hist.BinEdges(end)+mesh_pad;
    xlabel('|\omega|^2');
    ylabel('Count');
    title(['|\omega|^2 Spectrum & k-Means Centroids (wSteps = ' num2str(wSteps) ')']);
    hold on

    yl = ylim;
    for c = 1:nComponents
        plot([km_centroids(c), km_centroids(c)], yl,'Color',colorList(c,:),'LineWidth',2)
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
    %     t_nudge = 5;
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
    %     om_window_spec = repmat(om_sq, 1, wSteps);
    %         om_spec(:,stepSize*(k-1) + 1 : stepSize*(k-1) + wSteps) = om_window_spec;

        subplot(2,2,2);

        for g = 1:nComponents
            om_window_spec_cat = om_sq(om_class == g,:);
            if isempty(om_window_spec_cat) == 1
                continue
            end
            if logScale == 1
                for f = 1:length(om_window_spec_cat)
                    semilogy(mean(t),om_window_spec_cat(f),'.','Color',colorList(g,:),'LineWidth',1,'MarkerSize',10);
                end
            else
                for f = 1:length(om_window_spec_cat)
                    plot(mean(t),om_window_spec_cat,'.','Color',colorList(g,:),'LineWidth',1,'MarkerSize',10);
                end
            end
            hold on
        end
        title(['|\omega|^2 Spectra (r = ' num2str(r) ')']);
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
end