% clear variables; close all; clc

load('nBody_data.mat');
addpath('optdmd-master');
n = size(z,2)/4;

x = z(:,1:2*n);
% y = z(:,n+1:2*n);

figure
plot(t,x(:,1:n))
hold on
plot(t,x(:,n+1:2*n))
hold off
figure
plot(t,(x(:,1:n).^2 + x(:,n+1:2*n).^2).^(1/2))
title('r')
figure
plot(t,unwrap(angle(x(:,1:n)+sqrt(-1).*x(:,n+1:2*n))))
title('\theta')

x = x.';
TimeSpan = t;

r = 8; %rank to fit w/ optdmd
%  imode = 1, fit full data, slower
%  imode = 2, fit data projected onto first r POD modes
%      or columns of varargin{2} (should be at least r
%      columns in varargin{2})

nComponents = 2;

global_SVD = 1; %use global SVD modes for each DMD rather than individual SVD on each window

if global_SVD == 1
    [u,~,v] = svd(x,'econ');
    u = u(:,1:r); %just first r modes
    v = v(:,1:r);
end

initialize_artificially = 0; %this is required for use_last_freq OR use_median_freqs
use_last_freq = 0; 
use_median_freqs = 0; %mutually exclusive w/ use_last_freq

if use_median_freqs == 1
    load('km_centroids.mat');
    freq_meds = repelem(sqrt(km_centroids),r/nComponents);
end

nSplit = 10; %number of windows if they were non-overlapping

% stepSize = wSteps/10;
stepSize = 50;

thresh_pct = 1;

% save('mwDMD_params.mat','r','nComponents','initialize_artificially','use_last_freq','use_median_freqs','wSteps','nSplit','nSteps','nVars','thresh_pct','stepSize','nSlide');

%% execute optDMD

windows = [50 100 500 1000 1500];
wCounts = floor(length(TimeSpan)./windows);

for wc = 1:length(windows)
    wSteps = windows(wc);
    nSplit = wCounts(wc);
%     figure
    [mr_res, km_centroids] = run_optDMD(x,TimeSpan,r,wSteps,nSplit,stepSize,nComponents,thresh_pct,use_last_freq,use_median_freqs);
    plot_optDMD(mr_res,km_centroids,x,TimeSpan,wSteps,nSplit,stepSize,nComponents,thresh_pct);
end

return;

%% Visualize Modes
if exist('mr_res','var') == 0
    if nonlinear == 0
        try
            load('mwDMD_mr_res_i2_linear.mat');
        catch ME
            load('mwDMD_mr_res_linear.mat');
        end
    else
        try
            load('mwDMD_mr_res_i2.mat');
        catch ME
            load('mwDMD_mr_res.mat');
        end
    end
end

allModes = zeros(nSlide,r,r);
allFreqs = zeros(nSlide,r);
catList = flipud(perms(1:r));
modeCats = zeros(nSlide,r);
modeCats(1,:) = 1:r; %label using order of 1st modes
allModes(1,:,:) = mr_res{1}.w;
allFreqs(1,:) = mr_res{1}.Omega;
[~,catsAsc] = sort(mr_res{1}.Omega.*conj(mr_res{1}.Omega)); %category labels in ascending frequency order
for k = 2:nSlide
    w_old = squeeze(allModes(k-1,:,:));
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

figure('units','pixels','Position',[0 0 1366 768])
colorlist = {'k','r','b','g'};
w = squeeze(allModes(1,:,:));
wPlots = cell(r,r);
wTrails = cell(r,r);
trailLength = 1000; %in window steps
% Plot time series
subplot(2,4,5:8)
plot(TimeSpan(1:nSteps),x(:,1:nSteps),'LineWidth',1.5)
xlim([TimeSpan(1) TimeSpan(nSteps)])
hold on
lBound = plot([mr_res{1}.t(1) mr_res{1}.t(1)],ylim,'r-','LineWidth',2);
hold on
rBound = plot([mr_res{1}.t(end) mr_res{1}.t(end)],ylim,'r-','LineWidth',2);
hold on
% Plot 1st frame
for dim = 1:r
    subplot(2,4,dim)
    wi = w(dim,:);
    for j = 1:r
        wPlots{dim,j} = plot(real(wi(j)),imag(wi(j)),'o','Color',colorlist{j},'MarkerSize',7);
        hold on
        wTrails{dim,j} = plot(real(wi(j)),imag(wi(j)),'-','Color',colorlist{j},'LineWidth',0.1);
        hold on
        wTrails{dim,j}.Color(4) = 0.3; % 50% opacity
    end
    title(['Proj. Modes into Dimension ' num2str(dim)])
    axis equal
    xlim([-1 1])
    ylim([-1 1])
    xlabel('Real');
    ylabel('Imag');
    plot(xlim,[0 0],'k:')
    hold on
    plot([0 0],ylim,'k:')
    hold off
end
legend([wPlots{r,catsAsc(1)},wPlots{r,catsAsc(2)},wPlots{r,catsAsc(3)},wPlots{r,catsAsc(4)}],{'LF Mode 1','LF Mode 2','HF Mode 1','HF Mode 2'},'Position',[0.93 0.65 0.05 0.2])
%Plot subsequent frames
for k = 2:nSlide
    w = squeeze(allModes(k,:,:));
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

if use_median_freqs == 0
    save('mwDMD_allModes.mat','allModes','allFreqs');
else
    save('mwDMD_allModes_i2.mat','allModes','allFreqs');    
end

%% Times series of mode coordinates
figure('units','pixels','Position',[100 100 1366 768])
for j = 1:r
    subplot(r,2,2*j-1)
    plot(real(allModes(:,:,j)))
    title(['Re[w_' num2str(j) ']'])
%     legend('a','b','r','\theta')
    if j == r
        xlabel('Window #')
    end
    subplot(r,2,2*j)
    plot(imag(allModes(:,:,j)))
    title(['Im[w_' num2str(j) ']'])
%     legend('a','b','r','\theta')
    if j == r
        xlabel('Window #')
    end
end

%% Times series of modes in polar coordinates
figure('units','pixels','Position',[100 100 1366 768])
for j = 1:r
    cMode = squeeze(allModes(:,j,:));
    rMode = abs(cMode);
    thMode = unwrap(angle(cMode));
    subplot(r,2,2*j-1)
    plot(rMode)
    title(['Magnitudes of Modes in Dim. ' num2str(j)])
    legend({'LF','LF','HF','HF'})
%     legend('a','b','r','\theta')
    if j == r
        xlabel('Window #')
    end
    subplot(r,2,2*j)
    plot(thMode)
    title(['Angles of Modes in Dim. ' num2str(j)])
    legend({'LF','LF','HF','HF'})
%     legend('a','b','r','\theta')
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

%% Output mode time series data
% all_freqs_sq = allFreqs .* conj(allFreqs);
% med_high_band = median(max(all_freqs_sq,[],2));
% 
% ceil_pct = 1.05; %get rid of outliers w/ unnaturally high frequency
% 
% % ol_ind = max(all_freqs_sq,[],2) >=  ceil_pct * med_high_band;
% ol_ind = 628:642; %manually excise outliers from DMD error
% allFreqsCut = allFreqs;
% allFreqsCut(ol_ind,:) = NaN;
% 
% t_step = (TimeSpan(2)-TimeSpan(1))*stepSize; %time per window sliding step
% wSpan = 0:t_step:(nSlide-1)*t_step;
% 
% [F,TF] = fillmissing(allFreqsCut,'spline','SamplePoints',wSpan);
% 
% for j = 1:r
%     figure
%     plot(wSpan,real(allFreqsCut(:,j)),'.',wSpan(TF(:,j)),real(F(TF(:,j),j)),'o')
%     hold on
%     plot(wSpan,real(allFreqs(:,j)),'k:')
%     hold off
% %     xlim([6.1 6.6])
% end
% 
t_step = (TimeSpan(2)-TimeSpan(1))*stepSize; %time per window sliding step
start_steps = 100; %chop off from beginning to get rid of initial transients
good_winds = [start_steps:626, 644:nSlide]; %manually excise outliers from DMD error
cutoff_inds = 626-start_steps; %after omitting outliers, this is the index location of the cutoff between continuous regions

allFreqsCut = allFreqs(good_winds,:);
allModesCut = allModes(good_winds,:,:);

modeStack = squeeze(reshape(allModesCut,size(allModesCut,1),r*r,1));

if use_median_freqs == 0
    save('modeSeries.mat','modeStack','t_step','cutoff_inds');
else
    save('modeSeries_i2.mat','modeStack','t_step','cutoff_inds');
end

%% Split HF/LF Reconstruction
xr_H = zeros(size(x(:,1:nSteps)));
xr_L = zeros(size(x(:,1:nSteps)));
all_b = zeros(nVars,nSlide);
xn = zeros(nSteps,1); %count # of windows contributing to each step
for k = 1:nSlide
    w = mr_res{k}.w;
    b = mr_res{k}.b;
    all_b(:,k) = b;
    Omega = mr_res{k}.Omega;
    om_class = mr_res{k}.om_class;
    t = mr_res{k}.t;
    c = mr_res{k}.c;
    t_start = mr_res{k}.t_start;
    tShift = t-t(1); %compute each segment of xr starting at "t = 0"
%     t_nudge = 5;
    % constant shift gets put in LF recon
    xr_L_window = w(:, om_class == 1)*diag(b(om_class == 1))*exp(Omega(om_class == 1)*(t-t_start)) + c;
    xr_H_window = w(:, om_class == 2)*diag(b(om_class == 2))*exp(Omega(om_class == 2)*(t-t_start));
    xr_L(:,(k-1)*stepSize+1:(k-1)*stepSize+wSteps) = xr_L(:,(k-1)*stepSize+1:(k-1)*stepSize+wSteps) + xr_L_window;
    xr_H(:,(k-1)*stepSize+1:(k-1)*stepSize+wSteps) = xr_H(:,(k-1)*stepSize+1:(k-1)*stepSize+wSteps) + xr_H_window;
    xn((k-1)*stepSize+1:(k-1)*stepSize+wSteps) = xn((k-1)*stepSize+1:(k-1)*stepSize+wSteps) + 1;
end
xr_L = xr_L./repmat(xn.',nVars,1);
xr_H = xr_H./repmat(xn.',nVars,1);
figure
subplot(3,1,1)
plot(TimeSpan(1:nSteps),x(:,1:nSteps))
title('Input Data')
subplot(3,1,2)
plot(TimeSpan(1:nSteps),xr_L)
title('LF Recon')
subplot(3,1,3)
plot(TimeSpan(1:nSteps),xr_H)
title('HF Recon')
figure
plot(all_b.')

save('mwDMD_sep_recon.mat','xr_L','xr_H');

%% Construct A Matrices
if global_SVD ~= 1
    disp('Warning: Averaging A matrices in different bases')
end
A_avg = zeros(r);
A_proj_avg = zeros(r);
A_thresh = 100;
for k = 1:nSlide
    x_window = mr_res{k}.x;
    w = mr_res{k}.w;
    Omega = mr_res{k}.Omega;
    A = w*diag(Omega)*pinv(w);
    A_proj =  (u'*w)*diag(Omega)*(pinv(w)*u);
    
    A(abs(A) <= A_thresh) = 0; %threshold small values
    A_proj(abs(A_proj) <= A_thresh) = 0;
    
    A_avg = A_avg + A;
    A_proj_avg = A_proj_avg + A_proj;
end
A_avg = real(A_avg)/nSlide; %strips small imaginary residuals
A_proj_avg = real(A_proj_avg)/nSlide;

figure
plot(TimeSpan(1:nSteps),v(1:nSteps,:))
legend({'Mode 1','Mode 2','Mode 3','Mode 4'})
