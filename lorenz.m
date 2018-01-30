close all; clear all

tmax = 300;

sigma = 10;
beta = 8/3;
rho = 28;
gamma = 1;
epsilon = .2;
lambda = 3;


p = 1;
q = .5;

f = @(t,a) [-sigma*a(1) + sigma*a(2) + epsilon * abs(a(4)-lambda).^4; rho*a(1) - a(2) - a(1)*a(3); -beta*a(3) + a(1)*a(2); ...
    (a(4) - (a(4).^3)/3 - a(5) + gamma*(a(1).^2 + a(2).^2 + a(3).^2)^(1/2)/(a(1)).^(1/2)); a(4) + p - q*a(5)];
[t,a] = ode45(f,[0 tmax],[15.7,16.8,35.9 0 0]);     % Runge-Kutta 4th/5th order ODE solver
plot3(a(:,1),a(:,2),a(:,3))

tdiffs = [t(2:end)-t(1:end-1) t(2:end)-t(1:end-1) t(2:end)-t(1:end-1)];
a_t = (a(2:end,1:3) - a(1:end-1,1:3)) .* (tdiffs.^(-1));

figure
subplot(3,1,1)
plot(t,(a(:,1).^2 + a(:,2).^2 + a(:,3).^2).^(1/2))
title('Displacement')
subplot(3,1,2)
plot(t(1:end-1),(a_t(:,1).^2 + a_t(:,2).^2 + a_t(:,3).^2).^(1/2))
title('Speed')
subplot(3,1,3)
plot(t,epsilon*abs(a(:,4)-lambda).^4)
title('Lorenz Feedback From F-N')


figure
plot(t,gamma*(a(:,1).^2 + a(:,2).^2 + a(:,3).^2).^(1/2)./a(:,1).^(1/2))
title('FitzHugh-Nagumo Input')

figure
plot(t,a(:,4),t,a(:,5))
title('FitzHugh-Nagumo Variables')


return;


dt = .25;
stoptime = 500;

sigma = 10;
rho = 28;
beta = 8/3;

t = 0;
x = 0.01;
y = -0.02;
z = 0.03;

tic
nframes = 1;
figure
while t<stoptime
    xnew = sigma*(y-x);
    ynew = x*(rho-z) - y;
    znew = x*y - beta*z;
    x = xnew;
    y = ynew;
    z = znew;
    
    t = t+dt;
%     hi.CData = B;
%     ht.String = ['Time = ' num2str(t)];
%     drawnow limitrate
    scatter3(x,y,z)
    drawnow limitrate
    nframes = nframes+1;
end