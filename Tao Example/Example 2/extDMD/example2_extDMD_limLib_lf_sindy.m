clear; close all; clc

addpath('../SINDy_utils');

load('lf_sindy_data_2_ext_limLib.mat');
load('res_list_2_ext_limLib.mat')
sepTrial = 4;

x = lf_bt(:,isnan(lf_bt(1,:))==0);
nVars = size(x,1);
orig_norms = zeros(nVars,1);
for j = 1:nVars
    orig_norms(j) = norm(x(j,:));
    x(j,:) = x(j,:)/norm(x(j,:)); %normalize so b(t) and db/dt have equal magnitudes
end
TimeSpan = t_full(isnan(lf_bt(1,:))==0).';

polyorder = 3;
usesine = 0;
n = size(lf_bt,1);

h = TimeSpan(2)-TimeSpan(1);
%% Get Window Split Times


windSize = res_list(sepTrial,4);
nSplitGood = nnz(goodSplits);
nSteps = length(TimeSpan);
wEndList = windSize:windSize:(nSplitGood-1)*windSize;
wEndList = fliplr(wEndList);

%% compute Derivative 
xfull = x;
TimeSpanFull = TimeSpan;
clear('TimeSpan')

x = x(:,3:end-3);
tspan = TimeSpanFull;
dx = zeros(size(x));
for j=1:size(x,2)
    jf = j + 3;
%     dx(:,j) = (1/(12*h)) * (-xfull(:,jf-2) + xfull(:,jf+2) - 8*xfull(:,jf+1) + 8*xfull(:,jf-1));
    dx(:,j) = (1/(2*h)) * (xfull(:,jf+1) - xfull(:,jf-1));
end
% dx = dx + eps*randn(size(dx));   % add noise

x = [zeros(nVars,2) x zeros(nVars,3)]; %re-add beginning & end steps for index sanity
dx = [zeros(nVars,2) dx zeros(nVars,3)];

for endStep = wEndList
    x(:,endStep-3:endStep+3) = [];
    dx(:,endStep-3:endStep+3) = [];
    tspan(endStep-3:endStep+3) = [];
end

%re-remove padding @ beginning and end
x(:,1:2) = [];
x(:,end-2:end) = [];
dx(:,1:2) = [];
dx(:,end-2:end) = [];
tspan(1:2) = [];
tspan(end-2:end) = [];


x = x.';
dx = dx.';

x0 = x(1,:);


%% pool Data  (i.e., build library of nonlinear time series)
Theta = poolData(x,n,polyorder,usesine);
m = size(Theta,2);

%% compute Sparse regression: sequential least squares
lambdas = 10.^(2.5 : 0.1 : 4.8);
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
lambda = lambdas(7);      % lambda is our sparsification knob.
Xi = sparsifyDynamics(Theta,dx,lambda,n)



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

coeffsUsed = Xi(Xi~=0);
stringLibUsed = stringLib(Xi~=0);
for j = 1:length(coeffsUsed)
    disp([num2str(coeffsUsed(j)) ' ' stringLibUsed{j}]);
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
