function xr = example2_tdDMD_plot(x,TimeSpan,wSteps,nSplit,r,mr_res,kmList,thresh_pct,nComponents,downsample,outFile) %if outFile is empty, no export
    %% Plot MultiRes Results
    addpath('../../altmany-export_fig-9ac0917');
    
    x = x(:,1:downsample:end);
    TimeSpan = TimeSpan(1:downsample:end);
    
    pn = 1;
    wSteps = floor(wSteps/downsample);
    nSteps = wSteps * nSplit;
    logScale = 0;
    
    % figure('units','pixels','Position',[0 0 1366 2*768])

    % plotDims = [3 4]; %rows, columns of plot grid on screen at a given time
    plotDims = [1 4]; %rows, columns of plot grid on screen at a given time
    colorList = {'b','r','g','k','y'};
    
    x_PoT = x(:,1:nSteps);
    t_PoT = TimeSpan(1:nSteps);
    xr = zeros(size(x_PoT));
    %res_list: [pn, level #, nSplit, sampleSteps/nSplit]

    dupShift = 0.02; %multiplier to shift duplicate values so they can be visually distinguished

    nBins = 64;



    % for q = 1:size(res_list,1)
    % for q = [4 11 17]
    for q = 1
        figure('units','pixels','Position',[100 100 1200 400])
        om_spec = zeros(r,nSteps);
        omIm_spec = zeros(r,nSteps);
        b_spec = zeros(r,nSteps);
    %     scrollsubplot(plotDims(1),plotDims(2),[plotDims(2)*q-1, plotDims(2)*q]);
        subplot(plotDims(1),plotDims(2),[plotDims(2)-1, plotDims(2)]);
        plot(t_PoT,real(x_PoT),'k-') %plot ground truth
        xMax = max(max(abs(x_PoT)));

        C = kmList(q,:);

        ylim(1.5*[-xMax, xMax]);
        hold on

        all_om = [];

        for k = 1:nSplit
            try
                Omega = mr_res{pn,k}.Omega;
            catch ME
                continue
            end
            all_om = [all_om; Omega];
        end
        all_om_sq = all_om .* conj(all_om);

        all_om_sq = sort(all_om_sq);
        all_om_sq_thresh = all_om_sq(floor(thresh_pct*length(all_om_sq)));
        all_om_sq = all_om_sq(1:floor(thresh_pct*length(all_om_sq)));

        subplot(plotDims(1),plotDims(2),plotDims(2)-2);
        om_hist = histogram(all_om_sq,nBins);
        mesh_pad = 10;
%         bin_mesh = om_hist.BinEdges(1)-mesh_pad:0.5:om_hist.BinEdges(end)+mesh_pad;
        xlabel('|\omega|^2');
        ylabel('Count');
        hold on

        for g = 1:nComponents
            plot([C(g) C(g)],[0 max(om_hist.BinCounts)],'Color',colorList{g},'LineWidth',2)
            hold on
        end

        for k = 1:nSplit
            try
                w = mr_res{pn,k}.w;
            catch ME
                continue
            end
            b = mr_res{pn,k}.b;
            Omega = mr_res{pn,k}.Omega;
            c = mr_res{pn,k}.c;
            t_start = mr_res{pn,k}.t_start;

    %         if isempty(gmmList{q}) == 0
    %             om_class = mr_res{pn,j,k}.om_class;
    %         end
            try
                om_class = mr_res{pn,k}.om_class;
            catch ME
                continue
            end
            t = mr_res{pn,k}.t;
            t = t(1:downsample:end);
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
            om_window_spec = repmat(om_sq, 1, wSteps);
            om_spec(:,(k-1)*wSteps+1:k*wSteps) = om_window_spec;

            subplot(plotDims(1),plotDims(2),plotDims(2)-3);

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

            ylim([0 all_om_sq_thresh]);
    %         p_om_sq = semilogy(t,om_window_spec,'LineWidth',3);
    %         title(['Frequency Spectrum for ' num2str(wSteps) '-Step Window']);
            xlabel('t')
            xlim([t_PoT(1) t_PoT(end)]);
            ylabel('| \omega |^2')
            hold on

            xr_window = w*diag(b)*exp(Omega*(t-t_start)) + c;
            xr(:,(k-1)*wSteps+1:k*wSteps) = xr_window;


    %         scrollsubplot(plotDims(1),plotDims(2),plotDims(2)*q-1:plotDims(2)*q);
            subplot(plotDims(1),plotDims(2),plotDims(2)-1:plotDims(2));
            if rankDefFlag == 0
                p_xr = plot(t,real(xr_window),'LineWidth',2);
            else
                p_xr = plot(t,real(xr_window),'-.','LineWidth',1);
            end
            title([num2str(wSteps) '-Step Window (Frequencies ~' num2str(1/(wSteps*(t(2)-t(1)))) ' Hz)']);
            xlabel('t')
            xlim([t_PoT(1) t_PoT(end)]);
            ylabel('Re[x]')
            hold on


        end
        for k = 1:nSplit %plot dotted lines between time splits
            t = mr_res{pn,k}.t;
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
        try strlength(outFile)
%             export_fig outFile '-pdf';
            saveas(gcf,outFile)
            saveas(gcf,outFile,'png')
            close(gcf)
        catch ME
            close(gcf)
        end
    end
end