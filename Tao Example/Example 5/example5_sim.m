clc; clear; close all;

%% parameters

T=48;

x0 = [1 0];
tau1 = 15;
a = 0.7;
b = 0.8;
Iext = 0.5;

y0 = [0 1];
eta = 0; %dampling
epsilon = 1;
tau2 = 5;

%% RK4 integration of the mixed system

dt = 0.0001;

TimeSpan = 0:dt:T;   TimeSteps=length(TimeSpan)-1;

% x = zeros(4,TimeSteps+1);
% x(:,1) = [x1_0; x2_0; y1_0; y2_0];


[t1, x] = ode45(@(t1,x) rhs1(x,tau1,a,b,Iext),TimeSpan,x0);
[t2, y] = ode45(@(t2,y) rhs2(y,eta,epsilon,tau2),TimeSpan,y0);

tStep = mean(diff(t1))*4;
nSteps = ceil(T/tStep);
tN = 0:tStep:T;
tN = tN(1:nSteps); %match to nSteps

xOld = x;
yOld = y;
x = interp1(t1,xOld,tN); %enforce evenly spaced time steps
y = interp1(t2,yOld,tN); %enforce evenly spaced time steps
TimeSpan = tN;


%%

figure
subplot(2,1,1)
plot(TimeSpan, x);

title('x data');

xlabel('Time');

subplot(2,1,2)
plot(TimeSpan, y);
title('y data')
xlabel('Time');

%% 
uv = [x(:,1) y(:,1)];
save('raw_data_5_unmixed.mat','uv','TimeSpan')

%%
function dx = rhs1(x,tau,a,b,Iext)
    % FitzHugh-Nagumo Model
    v = x(1); w = x(2);
    vdot = v - (v^3)/3 - w + Iext;
    wdot = (1/tau) * (v + a - b*w);
    dx = [vdot; wdot];
end

function dy = rhs2(y,eta,epsilon,tau)
    % Unforced Duffing Oscillator
    p = y(1); q = y(2);
    pdot = q;
    qdot = (1/tau) * (-2*eta*q - p - epsilon*p^3);
    dy = [pdot; qdot];
end