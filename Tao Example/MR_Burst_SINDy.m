clear variables; close all; clc

addpath('../../Expanded SINDy')

% these may be overwritten for specific data cases
ODEorder = 1; % x'(t) = f(x) or x''(t) = f(x)
polyorder = 1:3;
usesine = 0;
differentiation_method = 1; % 1 = standard central-difference, 2 = TV Regularized
smoothing = 0; %toggle smoothing of input data
lowpass = 0;

%% Initialize: Example 3 (untransformed x,y variables)
% 
% % Use HF mrDMD reconstruction
% inFile = 'Example 3/mwDMD/mwDMD_sep_recon.mat';
% load(inFile);
% y = real(xr_H); %start with HF dynamics
% lambdas = 10.^(1 : 0.2 : 4);
% lambdaIdx = 6; %which lambda to integrate on (use positive or negative index)
% useSVD = 1;
% svdRank = 2;

% Use LF mrDMD reconstruction
inFile = 'Example 3/mwDMD/mwDMD_sep_recon.mat';
load(inFile);
y = real(xr_L); %start with HF dynamics
lambdas = 10.^(1.5 : 0.1 : 3.5);
lambdaIdx = [12 5]; %which lambda to integrate on (use positive or negative index)
useSVD = 1;
svdRank = 2;

% 
% % % Use original y variables
% % inFile = 'Example 3/raw_data_3_hiRes.mat';
% % load(inFile);
% % y = x(:,3:4).';
% % tspan = TimeSpan;
% % lambdas = 10.^(1 : 0.2 : 4);
% % lambdaIdx = 2; %which lambda to integrate on (use positive or negative index)
% % useSVD = 1;
% % svdRank = 2;
% % % 
% % 
load('Example 3/mwDMD/mwDMD_params.mat')
nTrunc = [wSteps wSteps]; %number of steps to chop off beginning/end to avoid transients and differentiation artifacts
burstSize = 10000; %points per burst sample
burstSep = 12000; %steps between starts of consecutive bursts


%% Initialize: Example 2 (transformed a, b, r,theta variables)
% % Use original a,b,r,theta variables
% inFile = 'Example 2/raw_data_2_hiRes.mat';
% load(inFile);
% y = x;
% tspan = TimeSpan;
% lambdas = 10.^(3 : 0.2 : 6);
% lambdaIdx = 6;
% useSVD = 1;
% svdRank = 4;

% % Use HF mrDMD reconstruction
% inFile = 'Example 2/mwDMD/mwDMD_sep_recon.mat';
% load(inFile);
% y = real(xr_H); %start with HF dynamics
% lambdas = 10.^(1 : 0.2 : 4);
% lambdaIdx = 12; %which lambda to integrate on (use positive or negative index)
% useSVD = 1;
% svdRank = 2;


% % Use LF mrDMD reconstruction
% inFile = 'Example 2/mwDMD/mwDMD_sep_recon.mat';
% load(inFile);
% y = real(xr_L); %start with HF dynamics
% lambdas = 10.^(2.5 : 0.1 : 4.5);
% lambdaIdx = [8 7]; %which lambda to integrate on (use positive or negative index)
% useSVD = 1;
% svdRank = 2;
% 
% % 
% load('Example 2/mwDMD/mwDMD_params.mat');
% nTrunc = [wSteps wSteps];
% burstSize = 10000; %points per burst sample
% burstSep = 12000; %steps between starts of consecutive bursts


% % Use DMD mode coordinates
% inFile = 'Example 2/modeSeries_i2.mat';
% load(inFile);
% which_vars = 1; %0 = all vars, 1 = just HF, 2 = just LF
% if which_vars == 1
%     modeStack = modeStack(:,1:8);
% elseif which_vars == 2
%     modeStack = modeStack(:,9:16);
% end
% y = real(modeStack.');
% tspan = t_step*(1:size(y,2));
% lambdas = 10.^(1 : 0.1 : 4);
% lambdaIdx = 12; %which lambda to integrate on (use positive or negative index)
% useSVD = 0;
% svdRank = 2;
% polyorder = 1;

% nTrunc = [0 0];
% burstSize = 100; %points per burst sample
% burstSep = 110; %steps between starts of consecutive bursts
% 

%% Initialize: 3 Body System Example
% inFile = '3 Body Example/mwDMD_sep_recon.mat';
% load(inFile);
% load('3 Body Example/mwDMD_params.mat');
% downSample = 100; %downsample factor; no need for superhigh resolution for the slow variables
% for j = 1:length(xr_sep)
%     xr_sep{j} = xr_sep{j}(:,1:downSample:end);
% end
% tspan = tspan(1:downSample:end);
% 
% nTrunc = ceil([wSteps wSteps]/downSample);
% burstSize = ceil(200000/downSample); %points per burst sample
% burstSep = ceil(250000/downSample); %steps between starts of consecutive bursts
% 
% % Orbital 
% y = real(xr_sep{1});
% % lambdas = 10.^(0.5 : 0.1 : 2.5);
% % lambdaIdx = [9 7 12 17]; %which lambda to integrate on (use positive or negative index)
% % lambdas = 10.^(0 : 0.1 : 2);
% % lambdaIdx = [13 9 8 4];
% lambdas = 10.^(-1 : 0.1 : 1);
% lambdaIdx = [17 16 11 13];
% useSVD = 1;
% svdRank = 4;
% lowpass = 1;

% % Saturn 
% y = real(xr_sep{2});
% lambdas = 10.^(1 : 0.2 : 3);
% lambdaIdx = 4; %which lambda to integrate on (use positive or negative index)
% useSVD = 1;
% svdRank = 2;

% % Jupiter 
% y = real(xr_sep{3});
% lambdas = 10.^(1 : 0.2 : 3);
% lambdaIdx = 7; %which lambda to integrate on (use positive or negative index)
% useSVD = 1;
% svdRank = 2;


%% Preprocess
y = y(:,nTrunc(1)+1:end-nTrunc(2));
y = y/max(max(abs(y)));
tspan = tspan(nTrunc(1)+1:end-nTrunc(2));

if lowpass == 1
    filtTrunc = ceil([50000 0]/downSample); %additional steps to cut off from either end after filtering
    fcutlow=0.0002;   %low cut frequency in Hz
    fs = 1/(tspan(2)-tspan(1));
    [b,a]=butter(3,(fcutlow/fs)*2,'low');

    yBand = zeros(size(y));
    for j = 1:size(y,1)
        yBand(j,:) = filter(b, a, y(j,:));
    end

    filtShift = ceil(17943/downSample);
    % % Determine optimal filter shift
    % shiftIterations = 5;
    % testShifts = 1:5;
    % for shiftIter = 1:shiftIterations
    %     bestNorm = norm(yBand-y);
    %     optJ = 0;
    %     for j = 1:length(testShifts)
    %         thisShift = round(10^testShifts(j));
    %         testNorm = norm(y(:,1:end-thisShift) - yBand(:,thisShift+1:end));
    % %         figure
    % %         plot(y(:,1:end-thisShift).','k')
    % %         hold on
    % %         plot(yBand(:,thisShift+1:end).','r')
    % %         hold off
    % %         disp([shiftIter, j])
    %         if testNorm < bestNorm
    %             optJ = j;
    %             bestNorm = testNorm;
    %         end
    %     end
    %     if shiftIter < shiftIterations
    %         testShiftStep = testShifts(2)-testShifts(1);
    %         testShifts = linspace(testShifts(optJ)-testShiftStep/2, testShifts(optJ)+testShiftStep/2, length(testShifts));
    %     end
    % end
    % filtShift = round(10^(testShifts(optJ)));
    % clear('testShifts','testNorm','bestNorm','shiftIter','shiftIterations','testShiftStep');

    yBand = yBand(:,filtShift+1:end);
    tspan = tspan(1:end-filtShift);
    % figure
    % plot(tspan,y(:,1:end-filtShift).','k')
    % hold on
    % plot(tspan,yBand.','r')
    % hold off
    y = yBand; clear('yBand');

    y = y(:,filtTrunc(1)+1:end-filtTrunc(2));
    tspan = tspan(filtTrunc(1)+1:end-filtTrunc(2));
end
    
if smoothing == 1
    % moving-average filter:
    figure
    plot(y.','k')
%     title('Pre-Filter')
    filterLength = 50000; %steps
    B = (1/filterLength)*ones(filterLength,1);
    for j = 1:size(y,1)
        y(j,:) = filter(B,1,y(j,:));
    end
    y = y(:,round(filterLength/2):end-(filterLength-round(filterLength/2)));
    tspan = tspan(round(filterLength/2):end-(filterLength-round(filterLength/2)));
    hold on
    plot(y.','r')
%     title('Post-Filter')
end

nBurst = floor(size(y,2)/burstSep);


if useSVD == 1
    [U,S,V] = svd(y,'econ');
    plot(cumsum(diag(S))/sum(diag(S)),'o-')
    hold on
    plot([svdRank svdRank],[0 1],'k:')
    hold off
    title('SVD Spectrum Cutoff')
    U = U(:,1:svdRank);
    S = S(1:svdRank,1:svdRank);
    V = V(:,1:svdRank);
%     x = V.';
    x = (S * V.');
    nVars = svdRank;
else
    x = y;
    nVars = size(x,1);
end



h = tspan(2)-tspan(1);

%% Normalize x
% % Normalize rows of x separately
% orig_norms = ones(nVars,1);
% for j = 1:nVars
%     orig_norms(j) = norm(x(j,:));
%     x(j,:) = x(j,:)/orig_norms(j); %normalize so rows of x have equal magnitudes
% end

% Normalize all rows of x together (i.e. preserve relative magnitudes)
orig_norms = repmat(norm(x),nVars,1);
for j = 1:nVars
    x(j,:) = x(j,:)/orig_norms(j);
end

% Don't normalize
% orig_norms = ones(nVars,1);

%% compute Derivative 
if differentiation_method == 1
    [tspan,x,dx] = xDeriv(tspan,x);
elseif differentiation_method == 2
    dx = zeros(size(x,2)+1,size(x,1));
    alpha = 1000;
    for j = 1:size(x,1)
        dx(:,j) = TVRegDiff(x(j,:),5,alpha,[],'small',[],h,1,1);
    end
    dx = dx(2:end,:);
    dx = dx.';
    % data, # iter, reg. parameter, start guess, scale, epsilon, timestep, plot, verbose
end

% if ODEorder == 2
%     [tspan,~,dx] = xDeriv(tspan,dx); %suspicious that this doesn't work properly
%     x = x(:,2:end-1);
% end

x = x.';
dx = dx.';
tspan = tspan.';

x0 = x(1,:);

% % Threshold x and dx for outliers
% nSigma = 4; %number of SDs away from mean to classify as outlier
% x_var = var(x);
% x_mean = mean(x);
% x_o = (x - repmat(x_mean,size(x,1),1)).^2 >= nSigma * repmat(x_var,size(x,1),1);
% x_omit = any((x - repmat(x_mean,size(x,1),1)).^2 >= nSigma * repmat(x_var,size(x,1),1),2);
% 
% dx_var = var(dx);
% dx_mean = mean(dx);
% dx_omit = any((dx - repmat(dx_mean,size(dx,1),1)).^2 >= nSigma * repmat(dx_var,size(dx,1),1),2);
% 
% omit_index = any([x_omit dx_omit],2);
% 
figure
subplot(2,1,1)
plot(tspan,real(x));
% hold on
% plot(tspan(omit_index),zeros(nnz(omit_index),1),'r*')
% hold off
title('Input Data')
subplot(2,1,2)
plot(tspan,real(dx));
% hold on
% plot(tspan(omit_index),zeros(nnz(omit_index),1),'r*')
% hold off
title(['Input Derivative (Order ' num2str(ODEorder) ')'])
% 
% x = x(~omit_index,:);
% dx = dx(~omit_index,:);
% tspan = tspan(~omit_index);

%% Optionally: Just train on first half of data; use second half for prediction comparison
halfTrain = 1;
if halfTrain == 1
    halfSteps = floor(size(x,1)/2);
    tspan2 = tspan(halfSteps+1:end);
    tspan = tspan(1:halfSteps);
    x2 = x(halfSteps+1:end,:);
    x = x(1:halfSteps,:);
    dx2 = dx(halfSteps+1:end,:);
    dx = dx(1:halfSteps,:);
    nBurst = floor(nBurst/2);
end

%% pool Data  (i.e., build library of nonlinear time series)
Theta = poolData(x,nVars,polyorder,usesine);
m = size(Theta,2);

%% Normalize columns of Theta
meanNorm = 0;
ThetaNorms = ones(size(Theta,2),1);
for tc = 1:size(Theta,2)
    ThetaNorms(tc) = norm(Theta(:,tc),1); %set 1-norm to 1
    meanNorm = meanNorm + ThetaNorms(tc);
    Theta(:,tc) = Theta(:,tc)/ThetaNorms(tc);
end
meanNorm = meanNorm/size(Theta,2);

%% Tune lambda on just the first burst
coeff_cts = zeros(length(lambdas),nVars);
for lj = 1:length(lambdas)
    testLambda = lambdas(lj);
    Xi = sparsifyDynamics(Theta(1:burstSize,:),dx(1:burstSize,:),testLambda,nVars);
    for li = 1:nVars
        coeff_cts(lj,li) = nnz(Xi(:,li));
    end
end

figure
for li = 1:nVars
    semilogx(lambdas,coeff_cts(:,li),'*-','LineWidth',2,'DisplayName',['# Terms: x' num2str(li)])
    hold on
end
title('Tuning \lambda');
xlabel('\lambda');
ylabel('# Nonzero Coefficients');
legend
ylim([0 max(max(coeff_cts))+1])
grid on
hold on

lambda = zeros(size(lambdaIdx));
for j = 1:length(lambdaIdx)
    if lambdaIdx(j) > 0
        lambda(j) = lambdas(lambdaIdx(j));
        cct = coeff_cts(lambdaIdx(j));
    else
        lambda(j) = lambdas(end+lambdaIdx(j));
        cct = coeff_cts(end+lambdaIdx(j));
    end
end

%highlight chosen lambda(s)
if length(lambda) == 1
    plot([lambda lambda], [0 max(max(coeff_cts))+1],'k:','LineWidth',2,'DisplayName','Chosen \lambda');
else
    set(gca,'ColorOrderIndex',1)
    for j = 1:length(lambda)
        plot([lambda(j) lambda(j)], [0 max(max(coeff_cts))+1],':','LineWidth',2,'DisplayName',['\lambda_' num2str(j)]);
        hold on
    end
end
hold off


%% recompute regression for chosen lambda
Xi_thresh = nBurst / 4; %minimum # of bursts that a coefficient must appear in
burst_Xis = zeros(nBurst,size(Theta,2),nVars);
for jb = 1:nBurst
    dxBurst = dx((jb-1)*burstSep+1 : (jb-1)*burstSep+burstSize,:);
    ThetaBurst = Theta((jb-1)*burstSep+1 : (jb-1)*burstSep+burstSize,:);
    XiBurst = sparsifyDynamics(ThetaBurst,dxBurst,lambda,nVars);
    burst_Xis(jb,:,:) = XiBurst./repmat(ThetaNorms,1,nVars);
end
Xi_counts = squeeze(sum(burst_Xis ~= 0, 1));

figure
imagesc(Xi_counts)
colorbar
set(gca,'CLim',[0 nBurst])

%smaller library with only the entries voted on above:
Theta = poolData(x,nVars,polyorder,usesine); %reset Theta to unnormalized version
ThetaUsed = Theta(:,any((Xi_counts > Xi_thresh),2));
if size(ThetaUsed,2) == 0
    disp('No terms survived voting!')
end

% XiUsed = ThetaUsed\dx;  % Least-squares
XiUsed = sparsifyDynamics(ThetaUsed,dx,lambda,nVars); % Sparse thresholding

XiUsedFull = zeros(size(XiBurst));
XiInd = 1;
for jl = 1:size(Theta,2)
    if any((Xi_counts(jl,:) > Xi_thresh))
        XiUsedFull(jl,:) = XiUsed(XiInd,:);
        XiInd = XiInd + 1;
    end
end
Xi = XiUsedFull.*(Xi_counts > Xi_thresh);
clear('XiInd','XiUsedFull','XiUsed','jl');



% Xi(Xi_counts > Xi_thresh) = XiUsed(:);
% Xi(Xi_counts < Xi_thresh) = 0;

%     %% integrate true and identified systems
%     options = odeset('RelTol',1e-10,'AbsTol',1e-10*ones(1,n));
% 
%     [tB,xB]=ode45(@(t,x)sparseGalerkin(t,x,Xi,polyorder,usesine),tspan,x0,options);  % approximate
% 
%     sindy_res{wn}.t_recon = tB;
%     sindy_res{wn}.x_recon = xB;
    
    

%% FIGURES!!

if halfTrain == 0
    xA = x;
    tA = tspan;
else
    xA = [x; x2];
    tA = [tspan; tspan2];
end
% Xi = Xi;

stringLib = libStringsFixed(nVars,polyorder,usesine).';
stringLib = repmat(stringLib, nVars, 1);

for nd = 1:nVars
    disp(['\dot{x' num2str(nd) '} = '])
    coeffsUsed = Xi(Xi(:,nd)~=0,nd);
    stringLibUsed = stringLib(nd,Xi(:,nd)~=0);
    for j = 1:length(coeffsUsed)
        disp([num2str(coeffsUsed(j)) ' ' stringLibUsed{j}]);
    end
    disp(' ') %line break
end

% options = odeset('RelTol',1e-7,'AbsTol',1e-7*ones(1,n));
options = odeset('RelTol',1e-6);

[tB,xB]=ode45(@(t,xB)sparseGalerkin(t,xB,Xi,polyorder,usesine),tA,x0,options);  % approximate

if size(x,2) >= 3
    figure
    dtA = [0; diff(tA)];
    plot_xA = plot3(xA(:,1),xA(:,2),xA(:,3),'r','LineWidth',1.5);
    hold on
    dtB = [0; diff(tB)];
    plot_xB = plot3(xB(:,1),xB(:,2),xB(:,3),'k','LineWidth',1.5);
    hold off
    plot_xA.Color(4) = 0.3; % opacity
    plot_xB.Color(4) = 0.3; % opacity
    xlabel('x_1','FontSize',13)
    ylabel('x_2','FontSize',13)
    zlabel('x_3','FontSize',13)
    l1 = legend('True','Identified');
    title('Manifolds: True vs. Identified')
end
% figure
% plot(tA,xA(:,1),'r','LineWidth',1.2)
% hold on
% plot(tA,xA(:,2),'r','LineWidth',1.2)
% plot(tB(1:10:end),xB(1:10:end,1),'k','LineWidth',1.2)
% hold on
% plot(tB(1:10:end),xB(1:10:end,2),'k','LineWidth',1.2)
% xlabel('Time')
% ylabel('State, x_k')
% legend('True x_1','True x_2','Identified x_1','Identified x_2')
% title('Time Series: True vs. Identified')


%% Plot Time Series
% [test_t, test_x] = ode45(@test_fn,[tspan(1) tspan(end)],[x(1,1) x(1,2)]);
% obtained_eps = ((abs(coeffsUsed(2)) * orig_norms(2)) / (abs(coeffsUsed(1)) * orig_norms(1))).^(-2);
xA_rescale = xA .* repmat((orig_norms.^(-1)).', size(xA,1),1);
xB_rescale = xB .* repmat((orig_norms.^(-1)).', size(xB,1),1);
try
    figure
    plot_ylim = [min(min([xA_rescale xB_rescale])) max(max([xA_rescale xB_rescale]))];
    subplot(2,1,1)
    if halfTrain == 0
        plot(tA,xA_rescale,'LineWidth',1)
        ylim(plot_ylim);
    else
        plot(tA(1:halfSteps),xA_rescale(1:halfSteps,:),'-','LineWidth',1)
        hold on
        ax = gca;
        ax.ColorOrderIndex = 1;
        plot(tA(halfSteps+1:end),xA_rescale(halfSteps+1:end,:),'--','LineWidth',1)
        hold off
        ylim(plot_ylim);
    end
    title('Input Data (Ground Truth)')% \epsilon = 0.01)')
    subplot(2,1,2)
    if halfTrain == 0
        plot(tB,xB_rescale,'LineWidth',1)
        ylim(plot_ylim);
    else
        plot(tB(1:halfSteps),xB_rescale(1:halfSteps,:),'-','LineWidth',1)
        hold on
        ax = gca;
        ax.ColorOrderIndex = 1;
        plot(tB(halfSteps+1:end),xB_rescale(halfSteps+1:end,:),'--','LineWidth',1)
        hold off
        ylim(plot_ylim);
    end
    title('SINDy Result')% (Obtained \epsilon = ' num2str(obtained_eps) ')'])

    % plot(test_t,test_x)
    % hold on
    % plot(test_t,0.01*sin(c1*test_t),'k','LineWidth',1.5)
    % title('Test')

    % function dydt = test_fn(t,x)
    % %     c1 = 9.5;
    % %     c2 = -c1;
    %     c1 = 1;
    %     c2 = -0.01^(-1); %true value
    %     dydt = [c1 * x(2); c2 * x(1)];
    % end
catch ME
    disp(['Integration failed at t = ' num2str(tB(end))])
end
return;

%% Animate Results


figure('units','pixels','Position',[0 0 1366 768])
% first frame
wA1 = xA(1,:);
wA1 = reshape(wA1,r,rr);
wB1 = xB(1,:);
wB1 = reshape(wB1,r,rr);
wPlotsA = cell(r,rr);
wTrailsA = cell(r,rr);
wPlotsB = cell(r,rr);
wTrailsB = cell(r,rr);
trailLength = 1000; %in window steps

subplot(2,r,r+1:2*r)
p_tsA = plot(tA,real(xA),'k','LineWidth',1.5);
hold on
p_tsB = plot(tB,real(xB),'r','LineWidth',1.5);
hold off
xlim([tA(1) tA(end)])
legend([p_tsA(1) p_tsB(1)],'Actual','SINDy Recon.')
% hold on
% lBound = plot([mr_res{1}.t(1) mr_res{1}.t(1)],ylim,'r-','LineWidth',2);
% hold on
% rBound = plot([mr_res{1}.t(end) mr_res{1}.t(end)],ylim,'r-','LineWidth',2);
% hold on

for dim = 1:r
    subplot(2,r,dim)
    wiA = wA1(dim,:);
    wiB = wB1(dim,:);
    for j = 1:rr
        wPlotsA{dim,j} = plot(real(wiA(j)),imag(wiA(j)),'o','Color','k','MarkerSize',7);
        hold on
        wPlotsB{dim,j} = plot(real(wiB(j)),imag(wiB(j)),'o','Color','r','MarkerSize',7);
        hold on
        wTrailsA{dim,j} = plot(real(wiA(j)),imag(wiA(j)),'-','Color','k','LineWidth',0.1);
        hold on
        wTrailsB{dim,j} = plot(real(wiB(j)),imag(wiB(j)),'-','Color','r','LineWidth',0.1);
        hold on
        wTrailsA{dim,j}.Color(4) = 0.3; % 30% opacity
        wTrailsB{dim,j}.Color(4) = 0.3;
    end
    title(['Proj. Modes into Dimension ' num2str(dim)])
    axis equal
    xlim([-0.1 0.1])
    ylim([-0.1 0.1])
    xlabel('Real');
    ylabel('Imag');
    plot(xlim,[0 0],'k:')
    hold on
    plot([0 0],ylim,'k:')
    hold off
end
% legend([wPlots{r,1},wPlots{r,2},wPlots{r,3},wPlots{r,4}],{'LF Mode 1','LF Mode 2','HF Mode 1','HF Mode 2'},'Position',[0.93 0.65 0.05 0.2])

for k = 2:length(tA)
    wA = xA(k,:);
    wA = reshape(wA,r,rr);
    wB = xB(k,:);
    wB = reshape(wB,r,rr);

    for dim = 1:r
%         subplot(4,4,dim)
        wiA = wA(dim,:);
        wiB = wB(dim,:);
        for j = 1:rr
            wPlotsA{dim,j}.XData = real(wiA(j));
            wPlotsA{dim,j}.YData = imag(wiA(j));
            wPlotsB{dim,j}.XData = real(wiB(j));
            wPlotsB{dim,j}.YData = imag(wiB(j));
            if k > trailLength
                wTrailsA{dim,j}.XData = [wTrailsA{dim,j}.XData(2:end) real(wiA(j))];
                wTrailsA{dim,j}.YData = [wTrailsA{dim,j}.YData(2:end) imag(wiA(j))];
                wTrailsB{dim,j}.XData = [wTrailsB{dim,j}.XData(2:end) real(wiB(j))];
                wTrailsB{dim,j}.YData = [wTrailsB{dim,j}.YData(2:end) imag(wiB(j))];
            else
                wTrailsA{dim,j}.XData = [wTrailsA{dim,j}.XData real(wiA(j))];
                wTrailsA{dim,j}.YData = [wTrailsA{dim,j}.YData imag(wiA(j))];
                wTrailsB{dim,j}.XData = [wTrailsB{dim,j}.XData real(wiB(j))];
                wTrailsB{dim,j}.YData = [wTrailsB{dim,j}.YData imag(wiB(j))];
            end
        end
    end
%     lBound.XData = [mr_res{k}.t(1) mr_res{k}.t(1)];
%     rBound.XData = [mr_res{k}.t(end) mr_res{k}.t(end)];
    pause(0.05)
end