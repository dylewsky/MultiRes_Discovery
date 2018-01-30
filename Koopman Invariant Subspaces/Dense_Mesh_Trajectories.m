close all; clear variables; clc

mu = -.05;
lambda = -1;
A = [mu 0 0; 0 lambda -lambda ; 0 0 2* mu ]; % Koopman linear dynamics
% A = randn(2,2);
% A = [A; [0 0]];
% A = [A, [0; 0; 0]];

lMesh = 0.2;
[X1m, X2m] = meshgrid(-2:lMesh:2 , -1:lMesh:4);
nPts = size(X1m,1)*size(X1m,2);

X1mV = reshape(X1m,nPts,1);
X2mV = reshape(X2m,nPts,1);

tspan = 0:.01:10;
tspan = tspan(1:end-1);
allTraj1 = zeros((length(tspan)-1) * nPts, 2);
allTraj2 = zeros((length(tspan)-1) * nPts, 2);
figure
for j = 1:nPts
    y0 = [X1mV(j); X2mV(j); X1mV(j).^2];
    [t,y] = ode45 (@(t,y)A*y,tspan ,y0);
    allTraj1((length(tspan)-1)*(j-1) + 1 : (length(tspan)-1)*j,:) = y(1:end-1,1:2);
    allTraj2((length(tspan)-1)*(j-1) + 1 : (length(tspan)-1)*j,:) = y(2:end,1:2);
    y1F = parab_arc_length(y(:,1));
    subplot(1,2,1)
    plot(y(:,1),y(:,2))
    hold on
    subplot(1,2,2)
    plot(y1F,y(:,2))
    hold on
end
subplot(1,2,1)
xlim([-2,2]);
ylim([-1,4]);
title('Original Nonlinear Dynamics')
xlabel('x_1')
ylabel('x_2')
subplot(1,2,2)
xlim(parab_arc_length([-2,2]));
ylim([-1,4]);
title('Flattened Koopman Inv. Manifold')
xlabel('m_1')
ylabel('m_2')

