clear; close all; clc

addpath('../../../SINDy_utils');

load('modeSeries.mat');
r = 4;

which_vars = 0; %0 = all vars, 1 = just LF, 2 = just HF

if which_vars == 1
    modeStack = modeStack(:,1:r*r/2);
elseif which_vars == 2
    modeStack = modeStack(:,r*r/2+1:end);
end

x = modeStack.';
nVars = size(x,1);
nChunks = length(cutoff_inds) + 1;
orig_norms = zeros(nVars,1);
for j = 1:nVars
    orig_norms(j) = norm(x(j,:));
    x(j,:) = x(j,:)/norm(x(j,:)); %normalize so b(t) and db/dt have equal magnitudes
end

TimeSpan = 0:t_step:(size(x,2)-1)*t_step;

polyorder = 1;
usesine = 0;
n = nVars;

h = t_step;

%% compute Derivative 
xfull = x;
TimeSpanFull = TimeSpan;
TimeSpan = [];
xCrop = [];
dxCrop = [];
for ch = 1:nChunks
    if ch == 1
        xChunk = x(:,1:cutoff_inds(ch));
        tChunk = TimeSpanFull(:,1:cutoff_inds(ch));
    elseif ch == nChunks
        xChunk = x(:,cutoff_inds(end)+1:end);
        tChunk = TimeSpanFull(:,cutoff_inds(end)+1:end);
    else
        xChunk = x(:,cutoff_inds(ch-1)+1:cutoff_inds(ch));
        tChunk = TimeSpanFull(:,cutoff_inds(ch-1)+1:cutoff_inds(ch));
    end
    xChunkCrop = xChunk(:,3:end-2);
    tChunkCrop = tChunk(:,3:end-2);
    dxChunk = (1/(12*h)) * (xChunk(:,1:end-4) - xChunk(:,5:end) + 8*xChunk(:,2:end-3) - 8*xChunk(:,4:end-1));
    xCrop = [xCrop xChunkCrop];
    dxCrop = [dxCrop dxChunk];
    TimeSpan = [TimeSpan tChunkCrop];
end

x = xCrop.';
dx = dxCrop.';
tspan = TimeSpan.';

x0 = x(1,:);
% figure
% plot(real(x));
% figure
% plot(real(dx));

%% pool Data  (i.e., build library of nonlinear time series)
Theta = poolData(x,n,polyorder,usesine);
m = size(Theta,2);

%% compute Sparse regression: sequential least squares
lambdas = 10.^(-1 : 0.1 : 1.2);
coeff_cts = zeros(size(lambdas));
for lj = 1:length(lambdas)
    testLambda = lambdas(lj);
    Xi = sparsifyDynamics(Theta,dx,testLambda,n);
    coeff_cts(lj) = nnz(Xi);
end
figure
semilogx(lambdas,coeff_cts,'*','LineWidth',2)
title('Tuning the Sparse Thresholding Parameter');
xlabel('\lambda');
ylabel('# Nonzero Coefficients');
grid on
if which_vars == 0
    lambda = lambdas(15);
elseif which_vars == 1
    lambda = lambdas(5); 
elseif which_vars == 2
    lambda = lambdas(15); 
end
Xi = sparsifyDynamics(Theta,dx,lambda,n);
figure
spy(Xi)
hold on
if which_vars == 0
    plot([0 0],[0.5 1.5],'w-',[0 0],[1.5 5.5],'k-',[0 0],[5.5 9.5],'r-',[0 0],[9.5 13.5],'b-',[0 0],[13.5 17.5],'g-','LineWidth',5)
    plot([1.5 5.5]-1,[0 0],'k-',[5.5 9.5]-1,[0 0],'r-',[9.5 13.5]-1,[0 0],'b-',[13.5 17.5]-1,[0 0],'g-','LineWidth',5)
elseif which_vars == 1
    plot([0 0],[0.5 1.5],'w-',[0 0],[1.5 5.5],'k-',[0 0],[5.5 9.5],'r-','LineWidth',5)
    plot([1.5 5.5]-1,[0 0],'k-',[5.5 9.5]-1,[0 0],'r-','LineWidth',5)
elseif which_vars == 2
    plot([0 0],[0.5 1.5],'w-',[0 0],[1.5 5.5],'b-',[0 0],[5.5 9.5],'g-','LineWidth',5)
    plot([1.5 5.5]-1,[0 0],'b-',[5.5 9.5]-1,[0 0],'g-','LineWidth',5)
end
    
ylabel('Coefficients')
xlabel('Dotted Variable')
% set(gca,'YTickLabel',stringLib(:,1))
% stringLibList = stringLib{:,1};
% gca.YTickLabel = stringLib(:,1);



%% integrate true and identified systems
options = odeset('RelTol',1e-10,'AbsTol',1e-10*ones(1,n));

[tB,xB]=ode45(@(t,x)sparseGalerkin(t,x,Xi,polyorder,usesine),tspan,x0,options);  % approximate

%% FIGURES!!
tA = tspan;
xA = x;

figure
dtA = [0; diff(tA)];
plot(xA(:,1),xA(:,2),'r.','LineWidth',.1);
hold on
dtB = [0; diff(tB)];
plot(xB(:,1),xB(:,2),'k--','LineWidth',1.2);
xlabel('x_1','FontSize',13)
ylabel('x_2','FontSize',13)
l1 = legend('True','Identified');

figure
plot(tA,xA(:,1),'r.','LineWidth',.1)
hold on
plot(tA,xA(:,2),'b.','LineWidth',.1)
plot(tB(1:10:end),xB(1:10:end,1),'k--','LineWidth',1.2)
hold on
plot(tB(1:10:end),xB(1:10:end,2),'k--','LineWidth',1.2)
xlabel('Time')
ylabel('State, x_k')
legend('True x_1','True x_2','Identified')

stringLib = libStringsFixed(nVars,polyorder,usesine).';
stringLib = repmat(stringLib, 1, nVars);

for nd = 1:nVars
    disp(['\dot{x' num2str(nd) '} = '])
    coeffsUsed = Xi(Xi(:,nd)~=0,nd);
    stringLibUsed = stringLib(Xi(:,nd)~=0,nd);
    for j = 1:length(coeffsUsed)
        disp([num2str(coeffsUsed(j)) ' ' stringLibUsed{j}]);
    end
end

%%
% [test_t, test_x] = ode45(@test_fn,[tspan(1) tspan(end)],[x(1,1) x(1,2)]);
% obtained_eps = ((abs(coeffsUsed(2)) * orig_norms(2)) / (abs(coeffsUsed(1)) * orig_norms(1))).^(-2);
figure
xA_rescale = xA .* repmat((orig_norms.^(-1)).', size(xA,1),1);
xB_rescale = xB .* repmat((orig_norms.^(-1)).', size(xB,1),1);
subplot(2,1,1)
plot(tA,xA_rescale,'.','LineWidth',1)
title('Input Data (Ground Truth)')% \epsilon = 0.01)')
subplot(2,1,2)
plot(tB,xB_rescale,'LineWidth',2)
title(['SINDy Result'])% (Obtained \epsilon = ' num2str(obtained_eps) ')'])

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

%% Animate Results

if which_vars == 0
    rr = r;
    colorlist = {'k','r','b','g'};
else
    rr = r/2; %reduced rank for hf/lf
    if which_vars == 1
        colorlist = {'k','r'};
    elseif which_vars == 2
        colorlist = {'b','g'};
    end
end

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