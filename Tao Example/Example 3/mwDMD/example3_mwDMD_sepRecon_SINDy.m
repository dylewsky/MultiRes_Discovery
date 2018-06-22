clear variables; close all; clc

addpath('../../../../sparsedynamics/utils/')

inFile = 'mwDMD_sep_recon.mat';
% dataLabel = 'Lorenz';
% inFile = [dataLabel '_sindy_input.mat'];

load(inFile);
x = real(xr_H); %choose which regime to use
x = x/max(max(abs(x))); %set max abs value to 1
figure
subplot(2,1,1)
plot(tspan,xr_L)
title('LF Reconstruction')
subplot(2,1,2)
plot(tspan,xr_H)
title('HF Reconstruction')

% x = x(1:2,:);

nVars = size(x,1);

polyorder = 1:3;
usesine = 0;

nTrunc = [1000 1000]; %number of steps to chop off beginning/end to avoid transients and differentiation artifacts


lambdas = 10.^(1 : 0.2 : 4);
% lambdas = 10.^(6 : 0.2 : 10);
lambdaIdx = -11; %which lambda to integrate on (use positive or negative index)



%     TimeSpan = 0:t_step:(size(x,2)-1)*t_step;
x = x(:,nTrunc(1)+1:end-nTrunc(2));
TimeSpan = tspan(nTrunc(1)+1:end-nTrunc(2));

%     ODE_order = 2;
figure
plot(TimeSpan,x)
title('Input Data');


h = TimeSpan(2)-TimeSpan(1);

%% Normalize x
% % Normalize rows of x separately
% orig_norms = ones(nVars,1);
% for j = 1:nVars
%     orig_norms(j) = norm(x(:,j));
%     x(:,j) = x(:,j)/orig_norms(j); %normalize so rows of x have equal magnitudes
% end

% % Normalize all rows of x together (i.e. preserve relative magnitudes)
% orig_norms = repmat(norm(x),nVars,1);
% for j = 1:nVars
%     x(:,j) = x(:,j)/orig_norms(j);
% end

% Don't normalize
orig_norms = ones(nVars,1);

%% compute Derivative 
xfull = x;
TimeSpanFull = TimeSpan;

%     xCrop = x(:,5:end-4);
%     dxCrop = (1/(12*h)) * (-x(:,1:end-4) + x(:,5:end) - 8*x(:,2:end-3) + 8*x(:,4:end-1));
%     dxCrop = dxCrop(:,3:end-2);
%     tCrop = TimeSpan(5:end-4);

xCrop = x(:,2:end-1);
dxCrop = (1/(2*h)) * (x(:,3:end) - x(:,1:end-2));
tCrop = TimeSpan(2:end-1);

x = xCrop.';
dx = dxCrop.';
tspan = tCrop.';


x0 = x(1,:);

% Threshold x and dx for outliers
nSigma = 4; %number of SDs away from mean to classify as outlier
x_var = var(x);
x_mean = mean(x);
x_o = (x - repmat(x_mean,size(x,1),1)).^2 >= nSigma * repmat(x_var,size(x,1),1);
x_omit = any((x - repmat(x_mean,size(x,1),1)).^2 >= nSigma * repmat(x_var,size(x,1),1),2);

dx_var = var(dx);
dx_mean = mean(dx);
dx_omit = any((dx - repmat(dx_mean,size(dx,1),1)).^2 >= nSigma * repmat(dx_var,size(dx,1),1),2);

omit_index = any([x_omit dx_omit],2);

figure
subplot(2,1,1)
plot(tspan,real(x));
hold on
plot(tspan(omit_index),zeros(nnz(omit_index),1),'r*')
hold off
title('Input Data')
subplot(2,1,2)
plot(tspan,real(dx));
hold on
plot(tspan(omit_index),zeros(nnz(omit_index),1),'r*')
hold off
title('Input Derivative')

x = x(~omit_index,:);
dx = dx(~omit_index,:);
tspan = tspan(~omit_index);




%% pool Data  (i.e., build library of nonlinear time series)
Theta = poolData(x,nVars,polyorder,usesine);
m = size(Theta,2);

%% Normalize columns of Theta
meanNorm = 0;
ThetaNorms = ones(size(Theta,2),1);
for tc = 1:size(Theta,2)
    meanNorm = meanNorm + norm(Theta(:,tc));
    ThetaNorms(tc) = norm(Theta(:,tc));
    Theta(:,tc) = Theta(:,tc)/ThetaNorms(tc);
end
meanNorm = meanNorm/size(Theta,2);

%% compute Sparse regression: sequential least squares
    
coeff_cts = zeros(length(lambdas),nVars);
for lj = 1:length(lambdas)
    testLambda = lambdas(lj);
    Xi = sparsifyDynamics(Theta,dx,testLambda,nVars);
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


if lambdaIdx > 0
    lambda = lambdas(lambdaIdx);
    cct = coeff_cts(lambdaIdx);
else
    lambda = lambdas(end+lambdaIdx);
    cct = coeff_cts(end+lambdaIdx);
end

%highlight chosen lambda
plot([lambda lambda], [0 max(max(coeff_cts))+1],'k:','LineWidth',2,'DisplayName','Chosen \lambda');
hold off

%% recompute regression for chosen lambda

Xi = sparsifyDynamics(Theta,dx,lambda,nVars);
Xi = Xi./repmat(ThetaNorms,1,nVars); %undo Theta normalization for coefficients



%     %% integrate true and identified systems
%     options = odeset('RelTol',1e-10,'AbsTol',1e-10*ones(1,n));
% 
%     [tB,xB]=ode45(@(t,x)sparseGalerkin(t,x,Xi,polyorder,usesine),tspan,x0,options);  % approximate
% 
%     sindy_res{wn}.t_recon = tB;
%     sindy_res{wn}.x_recon = xB;
    
    

%% FIGURES!!

tA = tspan;
xA = x;
% Xi = Xi;

stringLib = libStringsFixed(nVars,polyorder,usesine).';
stringLib = repmat(stringLib, 1, nVars);

for nd = 1:nVars
    disp(['\dot{x' num2str(nd) '} = '])
    coeffsUsed = Xi(Xi(:,nd)~=0,nd);
    stringLibUsed = stringLib(Xi(:,nd)~=0,nd);
    for j = 1:length(coeffsUsed)
        disp([num2str(coeffsUsed(j)) ' ' stringLibUsed{j}]);
    end
    disp(' ') %line break
end

% options = odeset('RelTol',1e-7,'AbsTol',1e-7*ones(1,n));
options = odeset('RelTol',1e-6);

[tB,xB]=ode45(@(t,xB)sparseGalerkin(t,xB,Xi,polyorder,usesine),tspan,x0,options);  % approximate

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
figure
xA_rescale = xA .* repmat((orig_norms.^(-1)).', size(xA,1),1);
xB_rescale = xB .* repmat((orig_norms.^(-1)).', size(xB,1),1);
subplot(2,1,1)
plot(tA,xA_rescale,'LineWidth',1)
title('Input Data (Ground Truth)')% \epsilon = 0.01)')
subplot(2,1,2)
plot(tB,xB_rescale,'LineWidth',1)
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