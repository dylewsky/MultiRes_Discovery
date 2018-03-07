clear all; close all; clc

%% Initialize
animate = 1;

rng(1234)
G = 6.674e-20; %km^3/(kg * s^2)

m_Sun = 1.989e30; %kg
m_Planet = 1.898e27;

n_Ast = 30;
m_Ast_expMean = 15;
m_Ast_expSD = 1;
m_Ast = 10.^(m_Ast_expMean + m_Ast_expSD * randn(n_Ast,1));
m_Bodies = [m_Sun; m_Planet; m_Ast];

n_Bodies = n_Ast+2;

x0_Sun = [0 0];
x0_Planet = [7.786e8 0]; %km
v0_Planet = [0 15]; %km/s

r0_Ast_mean = 4.0391e08;
r0_Ast_SD = 4.9866e07;
r0_Ast = r0_Ast_mean + r0_Ast_SD*randn(n_Ast,1);

th0_Ast = 2*pi*rand(n_Ast,1);
phi0_Ast = mod((th0_Ast + (pi/2)) + (pi/16)*randn(n_Ast,1),2*pi);

x0_Ast = r0_Ast .* [cos(th0_Ast) sin(th0_Ast)];
v0_Ast_mean = 22;
v0_Ast_SD = 1;
v0_Ast = (v0_Ast_mean + v0_Ast_SD*randn(n_Ast,1)) .* [cos(phi0_Ast) sin(phi0_Ast)];

p0_Ast = sum(repmat(m_Ast,1,2).*v0_Ast);
p0_Planet = m_Planet * v0_Planet;
p0_Sun = -(p0_Ast + p0_Planet);
v0_Sun = p0_Sun/m_Sun;

x0 = [x0_Sun; x0_Planet; x0_Ast];
x0 = [x0(:,1); x0(:,2)];
v0 = [v0_Sun; v0_Planet; v0_Ast];
v0 = [v0(:,1); v0(:,2)];


z0 = [x0; v0];

%% Simulate
tspan = [0 1e10];

allStable = 0;

nIter = 0;

while allStable == 0
    nIter = nIter + 1;
    disp(['Simulating: Iteration ' num2str(nIter)])
    [t, z] = ode45(@(t,z)nBody_RHS(z,m_Bodies,G),tspan,z0);

    d_final = (z(end,1:n_Bodies).^2 + z(end,n_Bodies+1:2*n_Bodies).^2).^(1/2);

    lBound = 1e8;
    uBound = 5e9;
    unst_bodies = (d_final < lBound) | (d_final > uBound);
    unst_bodies(1:2) = [0 0]; %keep sun and planet
    
    if nnz(unst_bodies) == 0
        allStable = 1;
    end

    m_Bodies(unst_bodies) = [];
    x0(repmat(unst_bodies,1,2)) = [];
    v0(repmat(unst_bodies,1,2)) = [];

    
    n_Bodies = n_Bodies - nnz(unst_bodies);
    
    %reset net momentum to 0
    p0_rocks = sum(repmat(m_Bodies(2:end),1,2).*[v0(2:n_Bodies) v0(n_Bodies+2:end)]);
    v0(1) = - p0_rocks(1)/m_Sun;
    v0(n_Bodies+1) = - p0_rocks(2)/m_Sun;
    
    z0 = [x0; v0];
end

%% Export Data
t_orig = t;
tStep = (tspan(2)-tspan(1))/size(z,1);
nSteps = ceil(tspan(2)/tStep);
tN = 0:tStep:tspan(2);
tN = tN(1:end-1); %match to nSteps
z_uneq = z;
z = interp1(t,z_uneq,tN);

% figure
% plot(t,r_uneq(:,1:3),'-o');
% hold on
% plot(tN,r(:,1:3),':.');
% hold off

t = tN;
save('nBody_data.mat','t','z');
% t = t_orig;
clear('t_orig','tN','tStep','nSteps');

%% Check Conservation Laws
x = z(:,1:n_Bodies);
y = z(:,n_Bodies+1:2*n_Bodies);
vx = z(:,2*n_Bodies+1:3*n_Bodies);
vy = z(:,3*n_Bodies+1:4*n_Bodies);
px = repmat(m_Bodies',length(t),1).*vx;
py = repmat(m_Bodies',length(t),1).*vy;
p_net = [sum(px,2) sum(py,2)];
figure
% plot(t,(p_net(:,1).^2 + p_net(:,2).^2).^(1/2))
plot(t,(p_net(:,1).^2 + p_net(:,2).^2).^(1/2)/sum(m_Bodies))
xlabel('t');
ylabel('COM Velocity (km/s)');

T = sum((1/2) * repmat(m_Bodies',length(t),1) .* (vx.^2 + vy.^2),2);
U = zeros(size(T));
for j = 1:n_Bodies
    for k = 1:j
        if j == k
            continue
        end
        dist = ((x(:,j) - x(:,k)).^2 + (y(:,j) - z(:,k)).^2).^(1/2);
        U = U - G*m_Bodies(j)*m_Bodies(k)./dist;
    end
end
figure
plot(t,T,t,U,t,T+U);
legend({'Kinetic', 'Potential', 'Total Energy'});


%% Plot Orbits
figure
sAlpha = 0.4;
x = z(:,1:n_Bodies);
y = z(:,n_Bodies+1:2*n_Bodies);

if animate == 0
    plot(x(:,1),y(:,1),'LineWidth',5); %sun
    axis equal
    hold on
    plot(x(:,2),y(:,2),'LineWidth',3); %planet
    hold on
    plot(x(:,3:end),y(:,3:end)) %asteroids
    hold off
    legend({'Sun','Planet','Asteroids'})
else
    sPlot = plot(x(1,1),y(1,1),'yo','MarkerSize',20,'MarkerFaceColor','y'); %sun
    hold on
    sTail = plot(x(1,1),y(1,1),'y-','LineWidth',10);
    sTail.Color = [1,1,0,sAlpha];
    axis equal
    hold on
    pPlot = plot(x(1,2),y(1,2),'ro','MarkerSize',10,'MarkerFaceColor','r'); %planet
    hold on
    pTail = plot(x(1,2),y(1,2),'r-','LineWidth',3);
    pTail.Color = [1,0,0,sAlpha];
    hold on
    for j = 1:n_Bodies - 2
        aPlot(j) = plot(x(1,3:end),y(1,3:end),'ko'); %asteroids
        hold on
        aTail(j) = plot(x(1,3:end),y(1,3:end),'k-');
        aTail(j).Color = [0,0,0,sAlpha];
    end
    
    hold off
    legend([sPlot pPlot aPlot(1)],{'Sun','Planet','Asteroids'})
    xlim([min(min(x)) max(max(x))]);
    ylim([min(min(y)) max(max(y))]);
    for k = 2:length(t)
        sPlot.XData = x(k,1);
        sPlot.YData = y(k,1);
        sTail.XData = x(1:k,1);
        sTail.YData = y(1:k,1);
        pPlot.XData = x(k,2);
        pPlot.YData = y(k,2);
        pTail.XData = x(1:k,2);
        pTail.YData = y(1:k,2);
        for j = 1:n_Bodies-2
            aPlot(j).XData = x(k,j);
            aPlot(j).YData = y(k,j);
            aTail(j).XData = x(1:k,j);
            aTail(j).YData = y(1:k,j);
        end
        pause(0.01);
    end
end

%% Plot Time Series
someAsteroids = randperm(n_Bodies-2,4)+2;
figure
x = z(:,1:n_Bodies);
y = z(:,n_Bodies+1:2*n_Bodies);
for j = 1:4
    subplot(2,2,j)
    plot(t,x(:,someAsteroids(j)));
    hold on
    plot(t,y(:,someAsteroids(j)));
    hold off
    xlim([t(1) t(end)]);
    title(['XY Pos, Body ' num2str(someAsteroids(j))])
end


