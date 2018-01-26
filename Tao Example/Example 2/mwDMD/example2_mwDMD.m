clear; close all; clc

% addpath('altmany-export_fig-9ac0917');
addpath('optdmd-master');
load('../raw_data_2_hiRes.mat');


r = size(x,1); %rank to fit w/ optdmd
imode = 1; %parameter for optdmd code
%  imode = 1, fit full data, slower
%  imode = 2, fit data projected onto first r POD modes
%      or columns of varargin{2} (should be at least r
%      columns in varargin{2})

nComponents = 2;
use_last_freq = 1;

wSteps = 11000;
nSplit = 20; %number of windows if they were non-overlapping
nSteps = wSteps * nSplit / 2;
nVars = size(x,1);
thresh_pct = 1;

% stepSize = wSteps/10;
stepSize = 100;

nSlide = floor((nSteps-wSteps)/stepSize);


%% execute optDMD

corner_sharpness = 64; %higher = sharper corners
lv_kern = tanh(corner_sharpness*(1:wSteps)/wSteps) - tanh(corner_sharpness*((1:wSteps)-wSteps)/wSteps) - 1;

mr_res = cell(nSlide);
for k = 1:nSlide
    sampleStart = stepSize*(k-1) + 1;
    sampleSteps = sampleStart : sampleStart + wSteps - 1;
    xSample = x(:,sampleSteps);
    tSample = TimeSpan(sampleSteps);

    mr_res{k}.x = xSample;
    xSample = xSample.*repmat(lv_kern,nVars,1); %round off corners
    mr_res{k}.t = tSample;
    
    c = mean(xSample,2);
    xSample = xSample - repmat(c,1,size(xSample,2));
    t_start = tSample(1);
    tSample = tSample - t_start;
    if (exist('e_init','var')) && (use_last_freq == 1)
        [w, e, b] = optdmd(xSample,tSample,r,imode,[],e_init);
    else
        [w, e, b] = optdmd(xSample,tSample,r,imode);
    end
    e_init = e;
    mr_res{k}.w = w;
    mr_res{k}.Omega = e;
    mr_res{k}.b = b;
    mr_res{k}.c = c;
    mr_res{k}.t_start = t_start;
end


%% Cluster Frequencies
close all;
if exist('mr_res','var') == 0
    load('mwDMD_mr_res.mat');
end

% nBins = 64;



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
save('mwDMD_mr_res.mat', 'mr_res');    


%% Plot MultiRes Results
close all;
if exist('mr_res','var') == 0
    load('mwDMD_mr_res.mat');
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
plot(t_PoT,real(x_PoT),'k-') %plot ground truth
xlabel('t')
ylabel('Measurement Data')
xMax = max(max(abs(x_PoT)));
ylim(1.5*[-xMax, xMax]);
hold on

xr = xr./repmat(xn.',4,1); %weight xr so all steps are on equal footing
plot(t_PoT,real(xr),'b-') %plot averaged reconstruction
title('Input (Black); DMD Recon. (Blue)');


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

%% Visualize Modes
if exist('mr_res','var') == 0
    load('mwDMD_mr_res.mat');
end

allModes = zeros(nSlide,r,r);
allFreqs = zeros(nSlide,r);
catList = flipud(perms(1:r));
modeCats = zeros(nSlide,r);
modeCats(1,:) = 1:r; %label using order of 1st modes
allModes(1,:,:) = mr_res{1}.w;
allFreqs(1,:) = mr_res{1}.Omega;
for k = 2:nSlide
    w_old = mr_res{k-1}.w;
    w = mr_res{k}.w;
    permDists = zeros(size(catList,1),1);
    for np = 1:size(catList,1)
        permDists(np) = norm(w(:,catList(np,:))-w_old);
    end
    [~,minCat] = min(permDists);
    modeCats(k,:) = catList(minCat,:);
    allModes(k,:,:) = w(:,catList(minCat,:));
    allFreqs(k,:) = mr_res{k}.Omega(catList(minCat,:));
end
meanFreqsSq = sum(allFreqs.*conj(allFreqs),1)/nSlide;

figure('units','pixels','Position',[100 100 1366 768])
colorlist = {'k','r','b','g'};
w = mr_res{k}.w;
wPlots = cell(r,r);
wTrails = cell(r,r);
trailLength = 10; %in steps
% Plot time series
subplot(4,4,5:16)
plot(TimeSpan(1:nSteps),x(:,1:nSteps),'LineWidth',1.5)
xlim([TimeSpan(1) TimeSpan(nSteps)])
hold on
lBound = plot([mr_res{1}.t(1) mr_res{1}.t(1)],ylim,'r-','LineWidth',2);
hold on
rBound = plot([mr_res{1}.t(end) mr_res{1}.t(end)],ylim,'r-','LineWidth',2);
hold on

% Plot 1st frame
for dim = 1:r
    subplot(4,4,dim)
    wi = w(dim,:);
    for j = 1:r
        wPlots{dim,j} = plot(real(wi(j)),imag(wi(j)),'o','Color',colorlist{j},'MarkerSize',7);
        hold on
        wTrails{dim,j} = plot(real(wi(j)),imag(wi(j)),'-','Color',colorlist{j},'LineWidth',0.2);
        hold on
    end
    title(['Proj. Modes into Dimension ' num2str(dim)])
    axis equal
    xlim([-1 1])
    ylim([-1 1])
    plot(xlim,[0 0],'k:')
    hold on
    plot([0 0],ylim,'k:')
    hold off
end
%Plot subsequent frames
for k = 2:nSlide
    w = mr_res{k}.w;
    for dim = 1:r
%         subplot(4,4,dim)
        wi = w(dim,:);
        for j = 1:r
            wPlots{dim,j}.XData = real(wi(j));
            wPlots{dim,j}.YData = imag(wi(j));
            if k > trailLength
                wTrails{dim,j}.XData = [wTrails{dim,j}.XData(2:end) real(wi(j))];
                wTrails{dim,j}.YData = [wTrails{dim,j}.YData(2:end) imag(wi(j))];
            else
                wTrails{dim,j}.XData = [wTrails{dim,j}.XData real(wi(j))];
                wTrails{dim,j}.YData = [wTrails{dim,j}.YData imag(wi(j))];
            end
        end
    end
    lBound.XData = [mr_res{k}.t(1) mr_res{k}.t(1)];
    rBound.XData = [mr_res{k}.t(end) mr_res{k}.t(end)];
    pause(0.05)
end

%% Times series of mode coordinates
figure('units','pixels','Position',[100 100 1366 768])
for j = 1:r
    subplot(r,2,2*j-1)
    plot(real(allModes(:,:,j)))
    title(['Re[w_' num2str(j) ']'])
    legend('a','b','r','\theta')
    if j == r
        xlabel('Window #')
    end
    subplot(r,2,2*j)
    plot(imag(allModes(:,:,j)))
    title(['Im[w_' num2str(j) ']'])
    legend('a','b','r','\theta')
    if j == r
        xlabel('Window #')
    end
end

%% Show how modes are reflections of one another
figure('units','pixels','Position',[100 100 1366 768])
subplot(2,2,1)
plot(real(allModes(:,:,1) - conj(allModes(:,:,2))))
title('Re[w_1 - w_2^*]')
subplot(2,2,2)
plot(imag(allModes(:,:,1) - conj(allModes(:,:,2))))
title('Im[w_1 - w_2^*]')

subplot(2,2,3)
plot(real(allModes(:,:,3) - conj(allModes(:,:,4))))
title('Re[w_3 - w_4^*]')
subplot(2,2,4)
plot(imag(allModes(:,:,3) - conj(allModes(:,:,4))))
title('Im[w_3 - w_4^*]')

