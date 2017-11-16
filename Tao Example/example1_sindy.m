clear; close all; clc

addpath('../SINDy_utils');

load('raw_data.mat');

polyorder = 5;
usesine = 1;
n = size(x,1);

h = TimeSpan(2)-TimeSpan(1);

%% compute Derivative 
xfull = x;
TimeSpanFull = TimeSpan;
clear('TimeSpan')
x = x(:,3:end-3);
tspan = TimeSpanFull(3:end-3);
dx = zeros(size(x));
for j=1:size(x,2)
    jf = j + 3;
    dx(:,j) = (1/(12*h)) * (xfull(:,jf-2) - xfull(:,jf+2) + 8*xfull(:,jf+1) - 8*xfull(:,jf-1));
%     dx(:,j) = (1/(2*h)) * (xfull(:,jf+1) - xfull(:,jf-1));
end
% dx = dx + eps*randn(size(dx));   % add noise

x = x.';
dx = dx.';

x0 = x(1,:);


%% pool Data  (i.e., build library of nonlinear time series)
Theta = poolData(x,n,polyorder,usesine);
m = size(Theta,2);

%% compute Sparse regression: sequential least squares
lambda = 0.05;      % lambda is our sparsification knob.
Xi = sparsifyDynamics(Theta,dx,lambda,n)

%% integrate true and identified systems
options = odeset('RelTol',1e-10,'AbsTol',1e-10*ones(1,n));

[tB,xB]=ode45(@(t,x)sparseGalerkin(t,x,Xi,polyorder,usesine),tspan,x0,options);  % approximate

%% FIGURES!!
tA = tspan;
xA = x;

figure
dtA = [0; diff(tA)];
plot(xA(:,1),xA(:,2),'r','LineWidth',1.5);
hold on
dtB = [0; diff(tB)];
plot(xB(:,1),xB(:,2),'k--','LineWidth',1.2);
xlabel('x_1','FontSize',13)
ylabel('x_2','FontSize',13)
l1 = legend('True','Identified');

figure
plot(tA,xA(:,1),'r','LineWidth',1.5)
hold on
plot(tA,xA(:,2),'b-','LineWidth',1.5)
plot(tB(1:10:end),xB(1:10:end,1),'k--','LineWidth',1.2)
hold on
plot(tB(1:10:end),xB(1:10:end,2),'k--','LineWidth',1.2)
xlabel('Time')
ylabel('State, x_k')
legend('True x_1','True x_2','Identified')

