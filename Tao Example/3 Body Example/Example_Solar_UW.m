%% Sun, Jupiter & Saturn dynamics
% Written for Daniel Dylewsky, Molei Tao and Nathan Kutz's project
% Jun 26 2018

clc;
clear;

%% initial conditions
global Gcst m0 m1 m2 epsilon

epsilon=0.001;              % Sun, Jupiter & Saturn
Gcst=4*pi^2;    m0=1;   m1=0.954791938424326609*epsilon;   m2=0.285885980666130812*epsilon;
daysPerYear=365.2422;
qIC=[[0;0;0], [4.84143144246472090E+00; -1.16032004402742839E+00; -1.03622044471123109E-01], [8.34336671824457987E+00; 4.12479856412430479E+00; -4.03523417114321381E-01]];
vIC=[[0;0;0], [1.66007664274403694E-03; 7.69901118419740425E-03; -6.90460016972063023E-05], [-2.76742510726862411E-03; 4.99852801234917238E-03;  2.30417297573763929E-05]]*daysPerYear;

paletteSize=10;     playMovie=0;    % parameters for animation

if (size(qIC,2)==3) % initial conditions include the star (star, planet1, planet2): convert to heliocentric coordiantes
    QIC=qIC;    PIC=[vIC(:,1)*m0 vIC(:,2)*m1 vIC(:,3)*m2];
    qIC=[QIC(:,2)-QIC(:,1), QIC(:,3)-QIC(:,1)];
    vIC=[PIC(:,2)/m1-(PIC(:,1)+PIC(:,2)+PIC(:,3))/(m0+m1+m2), PIC(:,3)/m2-(PIC(:,1)+PIC(:,2)+PIC(:,3))/(m0+m1+m2)];
end

T=100000;               % total simulation time, unit: year
h=0.1;  refinement=1;   % true timestep is h/refinement, but we only store every h for the sake of storage

%% Simulation in Cartesian by 4th-order symplectic integrator
tic

gamma=1/(2-2^(1/3));    % parameter for the 4th-order method
hh=h/refinement;

TimeSteps=floor(T/h+1e-8);
q=zeros(size(qIC,1),size(qIC,2),TimeSteps+1); q(:,:,1)=qIC;
v=zeros(size(qIC,1),size(qIC,2),TimeSteps+1); v(:,:,1)=vIC;
energy=zeros(1,TimeSteps+1);
i=0;
energy(i+1)=m1*norm(v(:,1,i+1))^2/2+m2*norm(v(:,2,i+1))^2/2-Gcst*m0*m1/norm(q(:,1,i+1))-Gcst*m0*m2/norm(q(:,2,i+1)) + norm(v(:,1,i+1)*m1+v(:,2,i+1)*m2)^2/2/m0 - Gcst*m1*m2/norm(q(:,1,i+1)-q(:,2,i+1));
for i=1:TimeSteps
    if mod(i*refinement,100000)==0
        disp([i*h T]);
    end
    
    qq=q(:,:,i);    vv=v(:,:,i);
    for refineIdx=1:refinement
        qq=qq+gamma*hh/2*(vv+(vv(:,1)*m1+vv(:,2)*m2)/m0*ones(1,2));
        vv=vv+gamma*hh*Gcst*(-[m0/norm(qq(:,1))^3*qq(:,1), m0/norm(qq(:,2))^3*qq(:,2)] - [m2/norm(qq(:,1)-qq(:,2))^3*(qq(:,1)-qq(:,2)), m1/norm(qq(:,1)-qq(:,2))^3*(qq(:,2)-qq(:,1))]);
        qq=qq+(1-gamma)/2*hh*(vv+(vv(:,1)*m1+vv(:,2)*m2)/m0*ones(1,2));
        vv=vv+(1-2*gamma)*hh*Gcst*(-[m0/norm(qq(:,1))^3*qq(:,1), m0/norm(qq(:,2))^3*qq(:,2)] - [m2/norm(qq(:,1)-qq(:,2))^3*(qq(:,1)-qq(:,2)), m1/norm(qq(:,1)-qq(:,2))^3*(qq(:,2)-qq(:,1))]);
        qq=qq+(1-gamma)/2*hh*(vv+(vv(:,1)*m1+vv(:,2)*m2)/m0*ones(1,2));
        vv=vv+gamma*hh*Gcst*(-[m0/norm(qq(:,1))^3*qq(:,1), m0/norm(qq(:,2))^3*qq(:,2)] - [m2/norm(qq(:,1)-qq(:,2))^3*(qq(:,1)-qq(:,2)), m1/norm(qq(:,1)-qq(:,2))^3*(qq(:,2)-qq(:,1))]);
        qq=qq+gamma*hh/2*(vv+(vv(:,1)*m1+vv(:,2)*m2)/m0*ones(1,2));
    end
    q(:,:,i+1)=qq;   v(:,:,i+1)=vv;

    energy(i+1)=m1*norm(v(:,1,i+1))^2/2+m2*norm(v(:,2,i+1))^2/2-Gcst*m0*m1/norm(q(:,1,i+1))-Gcst*m0*m2/norm(q(:,2,i+1)) + norm(v(:,1,i+1)*m1+v(:,2,i+1)*m2)^2/2/m0 - Gcst*m1*m2/norm(q(:,1,i+1)-q(:,2,i+1));
end

%% animate the result
disp(['Max energy fluctuation = ',num2str(max(abs(energy-energy(1))))]); % energy deviation (should be 0)
if (playMovie)
    figure
    surf([-paletteSize paletteSize], [-paletteSize paletteSize], [0 0; 0 0], 'EdgeColor','None','FaceAlpha',0.7,'FaceColor','green');
    hold on
    handle=scatter3(q(1,:,1),q(2,:,1),q(3,:,1),'filled');
    axis([-paletteSize paletteSize -paletteSize paletteSize -paletteSize/5 paletteSize/5]);
    coarse=1;
    for i=1:coarse:TimeSteps+1
        set(handle,'xdata',q(1,:,i),'ydata',q(2,:,i),'zdata',q(3,:,i));
        title(['Time=',num2str(h*(i-1))]);
        drawnow
    end
end

%% convert to orbital elements, which are slow variables
coarse=10;   % only convert every coarse data points (for increased speed)
aa=zeros(2,TimeSteps/coarse+1);
ee=zeros(2,TimeSteps/coarse+1);
ii=zeros(2,TimeSteps/coarse+1);
oomega=zeros(2,TimeSteps/coarse+1);
OOmega=zeros(2,TimeSteps/coarse+1);
EE=zeros(2,TimeSteps/coarse+1);
idx=0;
for i=1:coarse:TimeSteps+1
    if mod(i-1,10000*coarse)==0
        disp([i-1,TimeSteps]);
    end
    idx=idx+1;
    for j=1:2
        temp=computeOrbital(q(:,j,i),v(:,j,i));
        aa(j,idx)=temp(1);      ee(j,idx)=temp(2);
        ii(j,idx)=temp(3);      oomega(j,idx)=temp(4);
        OOmega(j,idx)=temp(5);  EE(j,idx)=temp(6);
    end
end

%% plot
% plot raw data (Cartesian coordinates)
figure
subplot(2,2,1);
plot([0:h:T], reshape(q(:,1,:),[3 size(q,3)]));
title('Jupiter position');
subplot(2,2,2);
plot([0:h:T], reshape(q(:,2,:),[3 size(q,3)]));
title('Saturn position');
subplot(2,2,3);
plot([0:h:T], reshape(v(:,1,:),[3 size(q,3)]));
title('Jupiter velocity');
subplot(2,2,4);
plot([0:h:T], reshape(v(:,2,:),[3 size(q,3)]));
title('Saturn velocity');


% plot underlying slow and fast data (orbital elements)
figure
subplot(2,3,1);
plot([0:h*coarse:T], aa);
title('a (slow)');
legend('Jupiter','Saturn');
subplot(2,3,2);
plot([0:h*coarse:T], ee);
title('e (slow)');
subplot(2,3,3);
plot([0:h*coarse:T], ii);
title('i (slow)');
subplot(2,3,4);
plot([0:h*coarse:T], unwrap(oomega(1,:)), [0:h*coarse:T], unwrap(oomega(2,:)));
title('\omega (slow)');
subplot(2,3,5);
plot([0:h*coarse:T], unwrap(OOmega(1,:)), [0:h*coarse:T], unwrap(OOmega(2,:)));
title('\Omega (slow)');
subplot(2,3,6);
plot([0:h*coarse:T], EE);
title('E (fast)');

%% Save Results
tspan = 0:h:T;
pos = [reshape(q(:,1,:),[3 size(q,3)]); reshape(q(:,2,:),[3 size(q,3)])];
vel = [reshape(v(:,1,:),[3 size(q,3)]); reshape(v(:,2,:),[3 size(q,3)])];
% Rows 1-3: Jupiter x,y,z
% Rows 4-6: Saturn x,y,z
save('Three_Body_Data_Cartesian.mat','tspan','pos','vel');

tcoarse = 0:h*coarse:T;
oomega = [unwrap(oomega(1,:)); unwrap(oomega(2,:))];
OOmega = [unwrap(OOmega(1,:)); unwrap(OOmega(2,:))];
save('Three_Body_Data_Orbital.mat','tcoarse','aa','ee','ii','oomega','OOmega','EE');
