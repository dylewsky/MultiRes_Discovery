clear; close all; clc

addpath('altmany-export_fig-9ac0917');
load('raw_data.mat');

nLevels = 5;
nVars = size(x,1);
nSteps = 2^14;
downScale = 3; %different-sized data sets will be built from blocks of size nSteps/2^downScale

if mod(nSteps * 2^(-downScale) * 2^(-(nLevels-1)),2) ~= 0
    print('Error: nSteps not sufficiently divisible by 2')
    return;
end

primeList = primes(2^downScale);
primeList = [primeList(2:end) 2^downScale]; %remove superfluous 2 and add 2^downScale

mr_res = cell(length(primeList),nLevels,2^(nLevels-1));
res_list = [];
for pn = 1:length(primeList)
    sampleSteps = nSteps * primeList(pn) / 2^(downScale);
    xSample = x(:,1:sampleSteps);
    tSample = TimeSpan(1:sampleSteps);
    
    nHold = 0;

    for n = 1:nLevels
        nSplit = 2^(n-1);
        xL = reshape(xSample, nVars, sampleSteps/nSplit, nSplit);
        tL = reshape(tSample, 1, sampleSteps/nSplit, nSplit);
        
        res_list = [res_list; pn, n, nSplit, sampleSteps/nSplit];
        for k = 1:nSplit
            xT = xL(:,:,k);
            tT = tL(:,:,k);
            mr_res{pn,n,k}.x = xT;
            mr_res{pn,n,k}.t = tT;
            [Phi, Omega, b] = run_dmd(xT,tT);
            mr_res{pn,n,k}.Phi = Phi;
            mr_res{pn,n,k}.Omega = Omega;
            mr_res{pn,n,k}.b = b;
        end
    %     nHold = input('Subtract off how many modes?')
    end
end

% mr_res_stack = reshape(mr_res,length(primeList)*nLevels*2^(nLevels-1),1,1);
res_list = sortrows(res_list,4,'descend');


%% Plot MultiRes Results
close all;

% figure('units','pixels','Position',[0 0 1366 2*768])

% plotDims = [3 4]; %rows, columns of plot grid on screen at a given time
plotDims = [1 4]; %rows, columns of plot grid on screen at a given time

x_PoT = x(:,1:nSteps);
t_PoT = TimeSpan(1:nSteps);
%res_list: [pn, level #, nSplit, sampleSteps/nSplit]

for q = 1:size(res_list,1)
    figure('units','pixels','Position',[100 100 1200 400])
    j = res_list(q,2);
    pn = res_list(q,1);
    nSplit = 2^(j-1);
    sampleSteps = nSteps * primeList(pn) / 2^(downScale);
    steps_per_window = sampleSteps/nSplit;
    om_spec = zeros(nVars,sampleSteps);
    b_spec = zeros(nVars,sampleSteps);

%     scrollsubplot(plotDims(1),plotDims(2),[plotDims(2)*q-1, plotDims(2)*q]);
    subplot(plotDims(1),plotDims(2),[plotDims(2)-1, plotDims(2)]);
    plot(t_PoT,real(x_PoT),'k-') %plot ground truth
    xMax = max(max(abs(x_PoT)));
    
    ylim(1.5*[-xMax, xMax]);
    hold on
    
    for k = 1:nSplit
        Omega = mr_res{pn,j,k}.Omega;
        Phi = mr_res{pn,j,k}.Phi;
        b = mr_res{pn,j,k}.b;
        t = mr_res{pn,j,k}.t;
        tShift = t-t(1); %compute each segment of xr starting at "t = 0"
        
        om_sq = conj(Omega).*Omega;
        om_window_spec = repmat(om_sq, 1, steps_per_window);
        om_spec(:,(k-1)*steps_per_window+1:k*steps_per_window) = om_window_spec;
        
%         scrollsubplot(plotDims(1),plotDims(2),plotDims(2)*q-3);
        subplot(plotDims(1),plotDims(2),plotDims(2)-3);
        plot(t,om_window_spec,'LineWidth',2);
%         title(['Frequency Spectrum for ' num2str(steps_per_window) '-Step Window']);
        xlabel('t')
        xlim([t_PoT(1) t_PoT(end)]);
        ylabel('| \omega |^2')
        hold on

%         if j ~=nLevels
%             set(gca,'XTick',[])
%         end
        
        b_sq = conj(b).*b;
        b_window_spec = repmat(b_sq, 1, steps_per_window);
%         scrollsubplot(plotDims(1),plotDims(2),plotDims(2)*q-2);
        subplot(plotDims(1),plotDims(2),plotDims(2)-2);
        plot(t,b_window_spec,'LineWidth',2);
%         title(['Weights for ' num2str(steps_per_window) '-Step Window']);
        xlabel('t')
        xlim([t_PoT(1) t_PoT(end)]);
        ylim auto;
        ylabel('|b|^2')
        hold on
        
%         if j ~=nLevels
%             set(gca,'XTick',[])
%         end
        
        xr_window = zeros(nVars,steps_per_window);
        for m = 1:nVars
            xr_window = xr_window + Phi(:,m) * exp(Omega(m)*tShift) * b(m);
        end
        
%         scrollsubplot(plotDims(1),plotDims(2),plotDims(2)*q-1:plotDims(2)*q);
        subplot(plotDims(1),plotDims(2),plotDims(2)-1:plotDims(2));
        plot(t,real(xr_window),'LineWidth',2);
        title([num2str(steps_per_window) '-Step Window (Frequencies ~' num2str(1/(steps_per_window*(t(2)-t(1)))) ' Hz)']);
        xlabel('t')
        xlim([t_PoT(1) t_PoT(end)]);
        ylabel('Re[x]')
        hold on

%         if j ~=nLevels
%             set(gca,'XTick',[])
%         end
    end
    for k = 1:nSplit %plot dotted lines between time splits
        t = mr_res{pn,j,k}.t;
%         scrollsubplot(plotDims(1),plotDims(2),plotDims(2)*q-3);
        subplot(plotDims(1),plotDims(2),plotDims(2)-3);
    	plot([t(end) t(end)],get(gca, 'YLim'),'k:')
        hold on 
%         scrollsubplot(plotDims(1),plotDims(2),plotDims(2)*q-2);
        subplot(plotDims(1),plotDims(2),plotDims(2)-2);
        plot([t(end) t(end)],get(gca, 'YLim'),'k:')
        hold on 
%         scrollsubplot(plotDims(1),plotDims(2),plotDims(2)*q-1:plotDims(2)*q);
        subplot(plotDims(1),plotDims(2),plotDims(2)-1:plotDims(2));
        plot([t(end) t(end)],get(gca, 'YLim'),'k:')
        hold on 
    end
    if q == 1
        export_fig 'manyres' '-pdf';
    else
        export_fig 'manyres' '-pdf' '-append';
    end
    close(gcf);
end

