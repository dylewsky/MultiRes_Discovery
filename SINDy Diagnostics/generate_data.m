clear variables; close all; clc

nStep = 990;
t_step = 0.01;

tspan = t_step*(0:nStep-1);

%% Slow oscillation

alpha = 10;
beta = 40;
delta = 0.02;

% x'' + delta x' + alpha * x + \beta x^3 = 0
% 
% x' = y
% y' = -delta* y - alpha * x - beta * x^3;

dxdt = @(x, alpha, beta, delta) [x(2); -delta * x(2) - alpha * x(1) - beta * x(1)^3];

x0 = rand(2,1);
[t1, y1] = ode45(@(t1,x)dxdt(x,alpha,beta,delta),tspan,x0);

plot(t1,y1(:,1))


%% Fast oscillation
omega = 22;

dydt = @(y, omega) [y(2); -omega^2 * y(1)];
y0 = rand(2,1);
[t2, y2] = ode45(@(t2,y)dydt(y,omega),tspan,y0);

plot(t2,y2(:,1))

%% Combined
if nnz(t1 ~= t2) > 0
    disp('Time vectors disagree')
end
epsilon = 20;

dzdt = @(z, alpha, beta, delta, omega) [z(2); -delta * z(2) - alpha * z(1) - beta * z(1)^3 - epsilon * z(3); z(4); -omega^2 * z(3)];

z0 = [x0; y0];

[t3, z] = ode45(@(t3,z)dzdt(z, alpha, beta, delta, omega),tspan,z0);

if nnz(t1 ~= t3) > 0
    disp('Time vectors disagree')
end

plot(t1,z(:,1))

%% Export data
t = t1;
save('testData.mat','y1','y2','z','t','alpha','beta','delta','omega','epsilon');

