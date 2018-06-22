clc; clear; close all;

%% parameters

a0=0;   b0=0;   r0=sqrt(2^2+0.8^2); theta0=atan2(0.8,2);

x1_0 = 0;
x2_0 = 0.5;
y1_0 = 0;
y2_0 = 0.5;

x0 = [x1_0; x2_0; y1_0; y2_0];

epsilon=0.01;
delta = 4;

T=48;



%% RK4 integration of the mixed system

h=epsilon/100;

TimeSpan = 0:h:T;   TimeSteps=length(TimeSpan)-1;

% x = zeros(4,TimeSteps+1);
% x(:,1) = [x1_0; x2_0; y1_0; y2_0];


[t, x] = ode45(@(t,x) rhs(x,delta,epsilon),TimeSpan,x0);

tStep = mean(diff(t))*4;
nSteps = ceil(T/tStep);
tN = 0:tStep:T;
tN = tN(1:nSteps); %match to nSteps

xOld = x;
x = interp1(t,xOld,tN); %enforce evenly spaced time steps
TimeSpan = tN;


%%

plot(TimeSpan, x);

title('Raw data');

xlabel('Time');


%% 
save('raw_data_3_hiRes.mat','x','TimeSpan')

%%
function dx = rhs(x,delta,epsilon)
    x1 = x(1); x2 = x(2); y1 = x(3); y2 = x(4);
    dx=[x2; -y1^2 * x1^3; y2; -epsilon^(-1)*y1 - delta*y1^3];
end